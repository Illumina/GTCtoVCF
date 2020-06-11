from vcf.parser import _Format
from math import log10

class GencallFormat(object):
    """
    Generate GQ format information for VCF
    """

    def __init__(self, logger, gencall_scores):
        """
        Create new GencallFormat object

        Args:
            logger (logging.Logger): Logging object for errors/warnings
            gencall_scores (list(float)): List of all gencall scores for all loci across sample

        Returns
            GencallFormat
        """
        self._logger = logger
        self._gencall_scores = gencall_scores

    @staticmethod
    def get_id():
        return "GQ"

    @staticmethod
    def get_description():
        return "GenCall score. For merged multi-locus entries, min(GenCall) score is reported."

    @staticmethod
    def get_format_obj():
        return _Format(
            GencallFormat.get_id(), 1, "Integer", GencallFormat.get_description())

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
            int: Minimum GenCall score, encoded as a phred quality (or None for empty arguments)
        """
        min_gencall_score = self._get_min_gencall_score(bpm_records)
        phred_quality_score = int(round(-10*log10(min(1.0, max(0.00001, 1 - min_gencall_score)))))
        return phred_quality_score
