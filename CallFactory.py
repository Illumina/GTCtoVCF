from collections import namedtuple
from vcf.model import _Call

class CallFactory(object):
    """
    The CallFactory is reponsible for generating _Calls to add to a VCF record
    """

    def __init__(self, formats, sample_name, logger):
        """
        Create new call factory

        Args:
            formats (.): List of Formats (GencallFormat, GenotypeFormat) which produce FORMAT values for a particular sample
            sample_name (string) : Name of sample reprsented by format
            logger (logging.logger): Logger for error logging

        Returns:
            CallFactory
        """
        self._format_namedtuple =  namedtuple("format_tuple", [format_element.get_id() for format_element in formats])
        self._formats = formats
        self._sample_name = sample_name
        self._logger = logger

    def create_call(self, locus_entry):
        """
        Update the VCF record with the call
        data from a new sample

        Args:
            locus_entry (LocusEntry) : Locus entry that links together BPM records to VCF record
            sample_name (string): Name of the sample

        Returns:
            vcf.model._Call - New call entry for this sample and locus
        """
        values = []

        for format_element in self._formats:
            values.append(format_element.generate_sample_format_info(
                locus_entry.bpm_records, locus_entry.vcf_record, self._sample_name))

        data = self._format_namedtuple(*values)
        return _Call(locus_entry.vcf_record, self._sample_name, data)

