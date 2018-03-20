from GenotypeFormat import GenotypeFormat
from GencallFormat import GencallFormat

class FormatFactory(object):
    """FormatFactory is responsible for determining the fields
    that will be present in the FORMAT string and creating
    the appropriate corresponding format objects
    """

    def __init__(self, no_samples, logger):
        """
        Create a new FormatFactory

        Args:
            no_samples (bool): True iff no sample data is present (format string should be empty)

        Returns:
            FormatFactory
        """
        self._logger = logger
        self._format_classes = []
        if not no_samples:
            self._format_classes.append(GenotypeFormat)
            self._format_classes.append(GencallFormat)

    def get_format_id_string(self):
        """
        Get the concatenated string of all format IDs

        Args:
            None

        Returns
            string: The concatenated string of all format IDs
        """
        return ":".join(self._get_format_ids())

    def _get_format_ids(self):
        """
        Get a list with all format IDs

        Args:
            None

        Returns:
            list(string): List of all format IDs
        """
        return [f.get_id() for f in self._format_classes]

    def get_format_objects(self):
        """
        Get a list with all formats (vcf.parser._Format)

        Args:
            None

        Returns
            list(vcf.parser._Format): List of all VCF formats
        """
        return [f.get_format_obj() for f in self._format_classes]

    def create_formats(self, gtc):
        """
        Create the list of formats

        Args:
            None

        Returns
            None
        """
        result = []
        result.append(GenotypeFormat(
            self._logger, gtc.get_gender(), gtc.get_genotypes()))
        result.append(GencallFormat(
            self._logger, gtc.get_genotype_scores()))
        return result
