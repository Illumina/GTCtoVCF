from vcf.parser import _Format

class GencallGenerator(object):
    """
    Generate GQ format information for VCF

    Attributes:
        id_string (string): "GQ"
        description (string): Format description
        format_obj (vcf._Format): VCF format object
    """

    def __init__(self, logger, gencall_scores):
        """
        Create new GenCallGenerator object

        Args:
            logger (logging.Logger): Logging object for errors/warnings
            gencall_scores (list(float)): List of all gencall scores for all loci across sample

        Returns
            GencallGenerator
        """
        self.id_string = "GQ"
        self.description = "GenCall score. For merged multi-locus entries, min(GenCall) score is reported."
        self.format_obj = _Format(
            self.id_string, 1, "String", self.description)
        self._logger = logger
        self._gencall_scores = gencall_scores

    def _get_min_gencall_score(self, bpm_records):
        """
        Look up gencall scores for entries defined by bpm_records and report
        the minimum.

        Args:
            bpm_records iter(BPMRecord) : Iterable of BPM records

        Returns
            float : Minimum gen call score
        """
        return min(self._gencall_scores[record.index_num] for record in bpm_records)

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """
        Returns sample format info (GQ score)

        Args:
            bpm_records (iter(BPMRecord)) : BPM records
            vcf_record (vcf._Record): VCF record
            sample_name (string): Name of sample

        Returns:
            float: Minimum GenCall score (or None for empty arguments)
        """
        return self._get_min_gencall_score(bpm_records)
