from .GenotypeFormat import GenotypeFormat
from .GencallFormat import GencallFormat
from .BAlleleFreqFormat import BAlleleFreqFormat
from .LogRRatioFormat import LogRRatioFormat


class FormatFactory(object):
    """FormatFactory is responsible for determining the fields
    that will be present in the FORMAT string and creating
    the appropriate corresponding format objects
    """

    def __init__(self, no_samples, formats_to_include, logger):
        """
        Create a new FormatFactory

        Args:
            no_samples (bool): True iff no sample data is present (format string should be empty)
            formats_to_include list(String): list of FORMAT IDs to include (default: GT,GC)

        Returns:
            FormatFactory
        """
        self._logger = logger
        self._format_classes = []
        for fmt in formats_to_include:
            assert fmt in FormatFactory.get_possible_formats(), "%r is not a valid FORMAT" % fmt
        self._formats_to_include = formats_to_include
        if not no_samples:
            if GenotypeFormat.get_id() in formats_to_include:
                self._format_classes.append(GenotypeFormat)
            if GencallFormat.get_id() in formats_to_include:
                self._format_classes.append(GencallFormat)
            if BAlleleFreqFormat.get_id() in formats_to_include:
                self._format_classes.append(BAlleleFreqFormat)
            if LogRRatioFormat.get_id() in formats_to_include:
                self._format_classes.append(LogRRatioFormat)

    @staticmethod
    def get_possible_formats():
        valid_classes = [GenotypeFormat.get_id(),
                         GencallFormat.get_id(),
                         BAlleleFreqFormat.get_id(),
                         LogRRatioFormat.get_id()]
        return valid_classes

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
        # only append if ID is in command line args
        if GenotypeFormat.get_id() in self._formats_to_include:
            result.append(GenotypeFormat(
                self._logger, gtc.get_gender(), gtc.get_genotypes()))
        if GencallFormat.get_id() in self._formats_to_include:
            result.append(GencallFormat(
                self._logger, gtc.get_genotype_scores()))
        if BAlleleFreqFormat.get_id() in self._formats_to_include:
            result.append(BAlleleFreqFormat(
                self._logger, gtc.get_ballele_freqs()))
        if LogRRatioFormat.get_id() in self._formats_to_include:
            result.append(LogRRatioFormat(
                self._logger, gtc.get_logr_ratios()))

        return result
