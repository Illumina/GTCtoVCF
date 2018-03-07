from operator import itemgetter
from LocusEntry import LocusEntry

class LocusEntryFactory(object):
    """
    Class to create locus entries from BPM records
    """

    def __init__(self, vcf_record_factory, locus_info_factory, call_factory, skip_indels, chrom_sort_function, unsquash_duplicates, logger):
        """
        Create new locus entry factory

        Args:
            vcf_record_factory (VcfRecordFactory): Creates vcf._Record objects
            locus_info_factory (LocusInfoFactory): Creates LocusInfo objects
            call_factory (CallFactory): Populates vcf._Record with _Call objects
            skip_indels (bool): True to skip indels
            chrom_sort_function (func(string, int)): Function used to sort chromosomes
            unsquash_duplicates (bool): True to generate separate entries for duplicates
            logger (logging.Logger): Logger to report warnings/errors
        """
        self._vcf_record_factory = vcf_record_factory
        self._locus_info_factory = locus_info_factory
        self._call_factory = call_factory
        self._chrom_sort_function = chrom_sort_function
        self._unsquash_duplicates = unsquash_duplicates
        self._logger = logger
        self._skip_indels = skip_indels

    def create_locus_entries(self, bpm_reader, loci_to_filter):
        """
        Generate locus entries from BPM records without
        any sample information.

        Args:
            bpm_reader (BPMReader): Provides BPM records
            loci_to_filter (set(string)): Set of loci names to filter from the manifest
        
        Returns:
            list(LocusEntry): List of locus entries corresponding to BPM records
        """
        result = []
        for record_group in self._group_bpm_records(bpm_reader.get_bpm_records(), loci_to_filter):
            result.append(self._generate_locus_entry(record_group))
        return sorted(result, key=lambda entry: (self._chrom_sort_function(entry.vcf_record.CHROM), entry.vcf_record.POS))

    def _group_bpm_records(self, bpm_records, loci_to_filter):
        """
        Group BPM records into groups where all BPM records in a single
        group will be represented in the same VCF record
        
        Args:
            bpm_reader (BPMReader): Provides BPM records
            loci_to_filter (set(string)): Set of loci names to filter from the manifest
        
        Yields:
            list(BPMRecord): Next group of BPM records
        """
        position2record = {}
        for record in bpm_records:
            if record.chromosome == "0" or record.pos == 0:
                continue
            if loci_to_filter is not None and record.name in loci_to_filter:
                continue
            if record.is_indel() and self._skip_indels:
                self._logger.warning("Skipping indel " + record.name)
                continue

            position = (record.chromosome, record.pos, None if record.indel_source_sequence is None else record.indel_source_sequence.get_plus_strand_sequence(record.ref_strand)[1])
            position2record.setdefault(position, []).append(record)
        for key, value in position2record.iteritems():
            if len(value) > 1 and self._unsquash_duplicates:
                alleles = set()
                for bpm_record in value:
                    alleles.update(bpm_record.get_plus_strand_alleles())
                if len(alleles) == 2:
                    for bpm_record in value:
                        yield [bpm_record]
                else:
                    yield value
            else:
                yield value

    def _generate_locus_entry(self, bpm_record_group):
        """
        Generate a single VCF record fromm a group of BPM records
        
        Args:
            bpm_record_group (list(BPMRecord)): Group of BPM records for single site
        
        Returns:
            LocusEntry: LocusEntry for the site
        """
        return LocusEntry(bpm_record_group, self._vcf_record_factory.create_vcf_record(bpm_record_group), self._locus_info_factory.create_locus_info(), self._call_factory)
