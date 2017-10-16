from InfoGenerator import LocusInfo

class LocusEntry(object):
    """
    Container class for BPM records and VCF record for single
    locus

    Attributes:
        bpm_records (list(BPMRecord)): List of BPM records for locus
        vcf_record (vcf._Record): VCF record for locus
    """

    def __init__(self, bpm_records, vcf_record, locus_info, call_factory):
        """
        Create new LocusEntry

        Args:
            bpm_records (list(BPMRecord)): Group of BPM records for a single locus
            vcf_record (vcf._Record): VCF record for a locus
            locus_info (LocusInfo): Info for the locus
            call_factory (CallFactory): Generates genotype information for locus

        Returns:
            LocusEntry
        """
        self.bpm_records = bpm_records
        self.vcf_record = vcf_record
        self._locus_info = locus_info
        self._call_factory = call_factory

    def add_sample(self, gtc, sample_name):
        """
        Update and both LocusInfo and call data from a sample GTC

        Args:
            gtc (GenotypeCalls): Genotyping information for sample
            sample_name (string): Sample name

        Returns:
            None
        """
        self._locus_info.add_sample(self.bpm_records, gtc)
        self._call_factory.add_sample_format(
            self.bpm_records, self.vcf_record, sample_name)

    def finalize(self):
        """
        Indicate no more samples are to be added and we can finalize the locus
        info in the VCF record

        Args:
            None

        Returns:
            None
        """
        self._locus_info.update(self.vcf_record)
