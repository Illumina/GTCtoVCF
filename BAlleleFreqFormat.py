from vcf.parser import _Format


class BAlleleFreqFormat(object):
    """
    Generate b allele frequency format information for VCF
    """
    def __init__(self, logger, b_allele_freq):
        self._b_allele_freq = b_allele_freq
        self._logger = logger

    @staticmethod
    def get_id():
        return "BAF"

    @staticmethod
    def get_description():
        return "B Allele Frequency"

    @staticmethod
    def get_format_obj():
        # arguments should be: ['id', 'num', 'type', 'desc']
        return _Format(BAlleleFreqFormat.get_id(), 1, "Float", BAlleleFreqFormat.get_description())

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """
        Get the sample b allele frequency

        Args:

        Returns:
            float: B allele frequency
        """
        # TODO: we need to handle the case where there are multiple BPMs
        idx = bpm_records[0].index_num
        return self._b_allele_freq[idx]
