from vcf.parser import _Format


class LogRRatioFormat(object):
    """
    Generate log R ratio format information for VCF
    """
    def __init__(self, logger, log_r_ratio):
        self._log_r_ratio = log_r_ratio
        self._logger = logger

    @staticmethod
    def get_id():
        return "LRR"

    @staticmethod
    def get_description():
        return "Log R Ratio"

    @staticmethod
    def get_format_obj():
        # arguments should be: ['id', 'num', 'type', 'desc']
        return _Format(LogRRatioFormat.get_id(), 1, "Float", LogRRatioFormat.get_description())

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """
        Get the sample log R ratio

        Args:

        Returns:
            float: log R ratio
        """
        # TODO: we need to handle the case where there are multiple BPMs
        idx = bpm_records[0].index_num
        return self._log_r_ratio[idx]
