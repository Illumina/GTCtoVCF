#!/usr/bin/env python2

import os
import sys
import logging
from argparse import ArgumentParser
import traceback
from vcf.parser import Writer, Reader
import tempfile

from BPMReader import BPMReader, CSVManifestReader
from InfoGenerator import LocusInfoFactory
from LocusEntryFactory import LocusEntryFactory
from CallFactory import CallFactory
from ReferenceGenome import ReferenceGenome
from IlluminaBeadArrayFiles import GenotypeCalls
from VcfRecordFactory import VcfRecordFactory
from ReaderTemplateFactory import ReaderTemplateFactory

VERSION = "1.0.1"

def is_dir_writable(parent_dir):
    try:
        test_file = tempfile.TemporaryFile(dir=parent_dir)
        test_file.close()
    except:
        return False
    return True

def is_file_writable(test_file):
    try:
        handle = open(test_file, "w")
        handle.close()
        os.remove(test_file)
    except:
        return False
    return True

def is_file_readable(test_file):
    try:
        handle = open(test_file, "r")
        handle.close()
    except:
        return False
    return True

def check_input_file(in_file, description, extensions):
    errors = []
    if not os.path.isfile(in_file):
        errors.append("Provided "+description+" file does not exist: "+in_file)
    else:
        if not is_file_readable(in_file):
            errors.append("Provided "+description+" is not readable: "+in_file)
        filename, extension = os.path.splitext(os.path.basename(in_file))
        if extensions and extension not in extensions:
            errors.append("Wrong file extension for "+description +" file provided: "+in_file +". Acceptable extensions are: "+",".join(extensions))
    return errors

def verify_inputs(args):
    """ Parse input arguments and report fatal errors
    """

    errors = []
    errors += check_input_file(args.manifest_file, 'manifest', ['.bpm', '.csv'])
    errors += check_input_file(args.genome_fasta_file, 'genome reference', [])

    if args.gtc_file:
        errors += check_input_file(args.gtc_file, 'gtc', ['.gtc'])

    is_curr_dir_writable = is_dir_writable(os.getcwd())

    if args.output_vcf_file:
        if not is_file_writable(args.output_vcf_file):
            errors.append("Cannot create file: " + args.output_vcf_file)
    else:
        if not is_curr_dir_writable:
            errors.append("Cannot write vcf to current directory: "+os.getcwd())

    return errors

def chrom_sort(chrom_string):
    """ Convert input chromosome string to appropriate type to ensure native
    	sort functions sorts chromosomes in order 1-22,X,Y,MT
    """
    return int(chrom_string) if str(chrom_string).isdigit() else chrom_string.lower()

def get_sample_name(gtc, gtc_file):

    sample_name = gtc.get_sample_name()
    if not sample_name:
        if ".gtc" in gtc_file:
            sample_name = (os.path.basename(gtc_file)).replace(".gtc", "")
        else:
            sample_name = "sample"
    return sample_name

def get_logger():
    # set up log file
    # create logger
    logger = logging.getLogger('GTC converter')
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def add_file_logger(logger, log_file):
    if log_file is None:
        return

    if not is_file_writable(log_file):
        raise Exception("Log file is not writable "+os.path.abspath(log_file))

    handler = logging.FileHandler(log_file)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def read_loci(loci_file):
    for line in open(loci_file):
        yield line.rstrip()

def read_auxiliary_records(auxiliary_loci):
    if auxiliary_loci is not None:
        auxiliary_records = {}
        with open(auxiliary_loci, "rb") as auxiliary_handle:
            for record in Reader(auxiliary_handle):
                auxiliary_records[record.ID] = record
        return auxiliary_records
    else:
        return None

def driver(gtc_file, manifest_file, genome_fasta_file, output_vcf_file, skip_indels, expand_identifiers, unsquash_duplicates, auxiliary_loci, loci_file, logger):
    loci_to_filter = set(read_loci(loci_file)) if loci_file is not None else None

    if manifest_file.lower().endswith(".bpm"):
        bpm_reader = BPMReader(manifest_file)
        if not skip_indels:
            raise Exception("Must skip indel processing (--skip-indels) when using BPM manifest, use CSV manifest to enable indel processing")
    elif manifest_file.lower().endswith(".csv"):
        bpm_reader = CSVManifestReader(manifest_file, logger)
    else:
        raise Exception("Manifest file must end with .bpm or .csv")

    sample_names = []
    if gtc_file:
        gtc = GenotypeCalls(gtc_file)
        if os.path.splitext(os.path.basename(gtc.get_snp_manifest()))[0].lower() != os.path.splitext(os.path.basename(manifest_file))[0].lower():
            raise Exception("Provided manifest name: "+manifest_file+ " and manifest file used to generate GTC file: "+gtc.get_snp_manifest()+" do not match. Exiting")
        logger.info("Manifest file used for GTC conversion identified as: " + gtc.get_snp_manifest())
        sample_name = get_sample_name(gtc, gtc_file)
        sample_names.append(sample_name)
    else:
        gtc = None
        logger.info("GTC file not provided")

    genome_reader = ReferenceGenome(genome_fasta_file, logger)

    locus_info_factory = LocusInfoFactory()

    call_factory = CallFactory(gtc, logger)
    
    reader_template = ReaderTemplateFactory(genome_reader, call_factory, "4.1", "gtc_to_vcf " + VERSION, chrom_sort, logger).create_reader_template(sample_names)

    auxiliary_records = read_auxiliary_records(auxiliary_loci)

    vcf_record_factory = VcfRecordFactory(call_factory, genome_reader, expand_identifiers, auxiliary_records, logger)

    locus_entries = LocusEntryFactory(vcf_record_factory, locus_info_factory, call_factory, skip_indels, chrom_sort, unsquash_duplicates, logger).create_locus_entries(bpm_reader, loci_to_filter)

    if gtc_file:
        for entry in locus_entries:
            if entry.vcf_record:
                entry.add_sample(gtc, sample_names[0])

    with open(output_vcf_file, "w") as output_handle:
        vcf_writer = Writer(output_handle, reader_template)
        for entry in locus_entries:
            if not entry.vcf_record:
                logger.warn("Could not create record for: "+entry.bpm_records[0].name)
                continue

            entry.finalize()
            vcf_writer.write_record(entry.vcf_record)

def main():
    parser = ArgumentParser(description="Convert GTC file to VCF format")
    parser.add_argument("--gtc-file", dest="gtc_file", required=False, help="GTC file")
    parser.add_argument("--manifest-file", dest="manifest_file", required=True, help="Bead pool manifest for product (*.csv or *.bpm)")
    parser.add_argument("--genome-fasta-file", dest="genome_fasta_file", required=True, help="Reference genome in fasta format")
    parser.add_argument("--output-vcf-file", dest="output_vcf_file", default="output.vcf", required=False, help="Name of VCF file to be created (default is output.vcf)")
    parser.add_argument("--skip-indels", dest="skip_indels", action="store_true", default=False, help="Skip processing of indels (default is False)")
    parser.add_argument("--log-file", dest="log_file", default=None, required=False, help="File to write logging information")
    parser.add_argument("--expand-identifiers", dest="expand_identifiers", action="store_true", default=False, help="For VCF entries with multiple corresponding manifest entries, list all manifest identifiers in VCF ID field")
    parser.add_argument("--unsquash-duplicates", dest="unsquash_duplicates", action="store_true", default=False, help="Generate unique VCF records for duplicate assays")
    parser.add_argument("--auxiliary-loci", dest="auxiliary_loci", default=None, required=False, help="VCF file with auxiliary definitions of loci")
    parser.add_argument("--filter-loci", dest="filter_loci", default=None, required=False, help="File containing list of loci names to filter from input manifest")
    parser.add_argument("--version", action="version", version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    logger = get_logger()
    add_file_logger(logger, args.log_file)
    errors = verify_inputs(args)

    if len(errors) > 0:
        for e in errors:
            logger.error(str(e))
        sys.exit(1)

    try:
        driver(args.gtc_file, args.manifest_file, args.genome_fasta_file, args.output_vcf_file, args.skip_indels, args.expand_identifiers, args.unsquash_duplicates, args.auxiliary_loci, args.filter_loci, logger)
    except Exception as e:
        logger.error(str(e))
        logger.debug(traceback.format_exc(e))
        sys.exit(-1)

if __name__ == "__main__":
    main()


