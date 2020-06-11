class LocusEntry(object):
    """
    Container class for BPM records and VCF record for single
    locus

    Attributes:
        bpm_records (list(BPMRecord)): List of BPM records for locus
        vcf_record (vcf._Record): VCF record for locus
    """

    def __init__(self, bpm_records, vcf_record):
        """
        Create new LocusEntry

        Args:
            bpm_records (list(BPMRecord)): Group of BPM records for a single locus
            vcf_record (vcf._Record): VCF record for a locus

        Returns:
            LocusEntry
        """
        self.bpm_records = bpm_records
        self.vcf_record = vcf_record

    def add_sample(self, call_factory, append):
        """
        Update and both LocusInfo and call data from a sample GTC

        Args:
            call_factory (CallFactory) : CallFactory to generate _Call objects
            append (bool) : If false, reset currrent calls beforing adding this sample

        Returns:
            None
        """
        call = call_factory.create_call(self)

        if self.vcf_record.samples is None or not append:
            self.vcf_record.samples = []
            self.vcf_record._sample_indexes = {}
        self.vcf_record._sample_indexes[call.sample] = len(self.vcf_record.samples)
        self.vcf_record.samples.append(call)
