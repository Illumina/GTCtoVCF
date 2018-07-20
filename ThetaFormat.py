from vcf.parser import _Format


class ThetaFormat(object):
    """
    Generate Theta format information for VCF
    """
    def __init__(self, logger, theta):
        self._theta = theta
        self._logger = logger

    @staticmethod
    def get_id():
        return "TH"

    @staticmethod
    def get_description():
        return "Normalized Theta-value for the sample"

    @staticmethod
    def get_format_obj():
        # arguments should be: ['id', 'num', 'type', 'desc']
        return _Format(ThetaFormat.get_id(), 1, "Float", ThetaFormat.get_description())

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """
        Get the sample genotype

        Args:

        Returns:
            float: Theta value
        """
        return self._theta
