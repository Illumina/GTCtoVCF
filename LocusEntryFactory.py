from LocusEntry import LocusEntry
from IlluminaBeadArrayFiles import RefStrand

class LocusEntryFactory(object):
    """
    Class to create locus entries from BPM records
    """

    def __init__(self, vcf_record_factory, chrom_sort_function, unsquash_duplicates, logger):
        """
        Create new locus entry factory

        Args:
            vcf_record_factory (VcfRecordFactory): Creates vcf._Record objects
            chrom_sort_function (func(string, int)): Function used to sort chromosomes
            unsquash_duplicates (bool): True to generate separate entries for duplicates
            logger (logging.Logger): Logger to report warnings/errors
        """
        self._vcf_record_factory = vcf_record_factory
        self._chrom_sort_function = chrom_sort_function
        self._unsquash_duplicates = unsquash_duplicates
        self._logger = logger

    def create_locus_entries(self, bpm_reader):
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
        for record_group in self._group_bpm_records(bpm_reader.get_bpm_records()):
            result.append(self._generate_locus_entry(record_group))
        return sorted(result, key=lambda entry: (self._chrom_sort_function(entry.vcf_record.CHROM), entry.vcf_record.POS))

    def _group_bpm_records(self, bpm_records):
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
            position = (record.chromosome, record.pos, (record.get_indel_source_sequences(RefStrand.Plus)[1], record.is_deletion) if record.is_indel() else None)
            position2record.setdefault(position, []).append(record)

        for _, value in position2record.iteritems():
            if len(value) > 1 and self._unsquash_duplicates:
                alleles = set()
                for bpm_record in value:
                    alleles.update(bpm_record.plus_strand_alleles)
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
        return LocusEntry(bpm_record_group, self._vcf_record_factory.create_vcf_record(bpm_record_group))
