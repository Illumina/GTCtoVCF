#!/usr/bin/env python2

import os
import sys
import logging
from argparse import ArgumentParser
import traceback
import tempfile
from vcf.parser import Writer, Reader

from BPMReader import BPMReader, CSVManifestReader
from LocusEntryFactory import LocusEntryFactory
from CallFactory import CallFactory
from ReferenceGenome import ReferenceGenome
from IlluminaBeadArrayFiles import GenotypeCalls
from VcfRecordFactory import VcfRecordFactory
from ReaderTemplateFactory import ReaderTemplateFactory
from FormatFactory import FormatFactory

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
        _, extension = os.path.splitext(os.path.basename(in_file))
        if extensions and extension not in extensions:
            errors.append("Wrong file extension for "+description +" file provided: "+in_file +". Acceptable extensions are: "+",".join(extensions))
    return errors

def verify_inputs(args):
    """ Parse input arguments and report fatal errors
    """
    errors = []
    errors += check_input_file(args.manifest_file, 'manifest', ['.bpm', '.csv'])
    errors += check_input_file(args.genome_fasta_file, 'genome reference', [])

    if args.gtc_files:
        for gtc_file in args.gtc_files:
            errors += check_input_file(gtc_file, 'gtc', ['.gtc'])

    if os.path.isdir(args.output_vcf_path):
        if not is_dir_writable(args.output_vcf_path):
            errors.append("Output directory " + args.output_vcf_path + " is not writeable")
    else:
        if not is_file_writable(args.output_vcf_path):
            errors.append("Unable to write to output file " + args.output_vcf_path)

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

def driver(gtc_files, manifest_file, genome_fasta_file, output_vcf_files, skip_indels, expand_identifiers, unsquash_duplicates, auxiliary_loci, loci_file, logger):
    loci_to_filter = set(read_loci(loci_file)) if loci_file is not None else None

    if manifest_file.lower().endswith(".bpm"):
        manifest_reader = BPMReader(manifest_file)
        if not skip_indels:
            raise Exception("Must skip indel processing (--skip-indels) when using BPM manifest, use CSV manifest to enable indel processing")
    elif manifest_file.lower().endswith(".csv"):
        manifest_reader = CSVManifestReader(manifest_file, logger)
    else:
        raise Exception("Manifest file must end with .bpm or .csv")

    auxiliary_records = read_auxiliary_records(auxiliary_loci)
    genome_reader = ReferenceGenome(genome_fasta_file, logger)
    format_factory = FormatFactory(gtc_files[0] is None, logger)
    reader_template_factory = ReaderTemplateFactory(genome_reader, format_factory, "4.1", "gtc_to_vcf " + VERSION, chrom_sort, logger)
    vcf_record_factory = VcfRecordFactory(format_factory, genome_reader, expand_identifiers, auxiliary_records, logger)
    locus_entries = LocusEntryFactory(vcf_record_factory, skip_indels, chrom_sort, unsquash_duplicates, logger).create_locus_entries(manifest_reader, loci_to_filter)

    for (gtc_file, output_vcf_file) in zip(gtc_files, output_vcf_files):
        if gtc_file:
            logger.info("Handling GTC file " + gtc_file)
            gtc = GenotypeCalls(gtc_file)
            if os.path.splitext(os.path.basename(gtc.get_snp_manifest()))[0].lower() != os.path.splitext(os.path.basename(manifest_file))[0].lower():
                logger.warn("Provided manifest name: "+manifest_file+ " and manifest file used to generate GTC file: "+gtc.get_snp_manifest()+" do not match, skipping")
                continue
            logger.info("Manifest file used for GTC conversion identified as: " + gtc.get_snp_manifest())
            sample_name = get_sample_name(gtc, gtc_file)
            reader_template = reader_template_factory.create_reader_template([sample_name])

            call_factory = CallFactory(format_factory.create_formats(gtc), sample_name, logger)
            for entry in locus_entries:
                if entry.vcf_record:
                    entry.add_sample(call_factory, False)
        else:
            logger.info("GTC file not provided")
            reader_template = reader_template_factory.create_reader_template([])

        with open(output_vcf_file, "w") as output_handle:
            vcf_writer = Writer(output_handle, reader_template)
            for entry in locus_entries:
                if not entry.vcf_record:
                    logger.warn("Could not create record for: "+entry.bpm_records[0].name)
                    continue
                vcf_writer.write_record(entry.vcf_record)

def generate_io_files(gtc_files, output_vcf_path, manifest_file):
    if not gtc_files:
        gtc_files = [None]
    
    if os.path.isdir(output_vcf_path):
        if gtc_files[0] is None:
            output_vcf_files = [os.path.join(output_vcf_path, os.path.splitext(os.path.basename(manifest_file))[0] + ".vcf")]
        else:
            output_vcf_files = []
            for gtc_file in gtc_files:
                output_vcf_files.append([os.path.join(output_vcf_path, os.path.splitext(os.path.basename(gtc_file))[0] + ".vcf"))
    else:
        if len(gtc_files) > 1:
            raise Exception("Must specify output directory (not file) when providing more than one input GTC file")
        output_vcf_files.append(output_vcf_path)

    return (gtc_files, output_vcf_files)


def main():
    parser = ArgumentParser(description="Convert GTC file to VCF format")
    parser.add_argument("--gtc-files", dest="gtc_files", nargs="+", required=False, help="GTC file")
    parser.add_argument("--manifest-file", dest="manifest_file", required=True, help="Bead pool manifest for product (*.csv or *.bpm)")
    parser.add_argument("--genome-fasta-file", dest="genome_fasta_file", required=True, help="Reference genome in fasta format")
    parser.add_argument("--output-vcf-path", dest="output_vcf_path", default="output.vcf", required=False, help="Name of VCF file to be created (default is output.vcf)")
    parser.add_argument("--skip-indels", dest="skip_indels", action="store_true", default=False, help="Skip processing of indels (default is False)")
    parser.add_argument("--log-file", dest="log_file", default=None, required=False, help="File to write logging information")
    parser.add_argument("--expand-identifiers", dest="expand_identifiers", action="store_true", default=False, help="For VCF entries with multiple corresponding manifest entries, list all manifest identifiers in VCF ID field")
    parser.add_argument("--unsquash-duplicates", dest="unsquash_duplicates", action="store_true", default=False, help="Generate unique VCF records for duplicate assays")
    parser.add_argument("--auxiliary-loci", dest="auxiliary_loci", default=None, required=False, help="VCF file with auxiliary definitions of loci")
    parser.add_argument("--filter-loci", dest="filter_loci", default=None, required=False, help="File containing list of loci names to filter from input manifest")
    parser.add_argument("--version", action="version", version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    args.output_vcf_path = os.path.abspath(args.output_vcf_path)
    args.manifest_file = os.path.abspath(args.manifest_file)
    args.genome_fasta_file = os.path.abspath(args.genome_fasta_file)

    logger = get_logger()
    add_file_logger(logger, args.log_file)
    errors = verify_inputs(args)

    if errors:
        for error in errors:
            logger.error(str(error))
        sys.exit(-1)

    (gtc_files, output_vcf_files) = generate_io_files(args.gtc_files, args.output_vcf_path, args.manifest_file)

    if gtc_files[0] is not None:
        logger.info("Processing %s GTC files", str(len(gtc_files)))
    try:
        driver(gtc_files, args.manifest_file, args.genome_fasta_file, output_vcf_files, args.skip_indels, args.expand_identifiers, args.unsquash_duplicates, args.auxiliary_loci, args.filter_loci, logger)
    except Exception as exception:
        logger.error(str(exception))
        logger.debug(traceback.format_exc(exception))
        sys.exit(-1)

if __name__ == "__main__":
    main()


