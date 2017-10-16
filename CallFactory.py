from collections import namedtuple
from vcf.model import _Call
from GenotypeGenerator import GenotypeGenerator
from GencallGenerator import GencallGenerator


class CallFactory(object):
    """
    The CallFactory is reponsible for adding and merging Calls to a VCF record
    """

    def __init__(self, gtc, logger):
        """
        Create new call factory

        Args:
            gtc (GenotypeCalls): Genotyping data for a sample
            logger (logging.logger): Logger for error logging

        Returns:
            CallFactory
        """
        self._gtc = gtc
        self._logger = logger

        if self._gtc:
            self._format_generators = self._create_format_generators()
            self._format_namedtuple = namedtuple(
                "format_tuple", self._get_format_ids())
        else:
            self._format_namedtuple = namedtuple("format_tuple", [])
            self._format_generators = []

    def format_cmp(self, format_string, _):
        """
        Comparsion function for FORMAT fields. Used to ensure
        that GT will be given first, as noted in the VCF specification.

        Args:
            format_string (string): Format string (e.g, "GT")

        Returns
            int: -1, 0, 1 for less than, equal, or greater
        """
        # If GT is available, put this at the top of the list.  Other elements
        # stay in same relative order
        if format_string.id_string == "GT":
            return -1
        return 1

    def get_format_id_string(self):
        """
        Get the concatenated string of all format IDs

        Args:
            None

        Returns
            string: The concatenated string of all format IDs
        """
        return ":".join([f.format_obj.id for f in self._format_generators])

    def _get_format_ids(self):
        """
        Get a list with all format IDs

        Args:
            None

        Returns:
            list(string): List of all format IDs
        """
        return [f.format_obj.id for f in self._format_generators]

    def get_format_objects(self):
        """
        Get a list with all format generators

        Args:
            None

        Returns
            list(Generator): List of all format generators
        """
        return [f.format_obj for f in self._format_generators]

    def _create_format_generators(self):
        """
        Create the list of format generators

        Args:
            None

        Returns
            None
        """
        format_generators = []
        format_generators.append(GencallGenerator(
            self._logger, self._gtc.get_genotype_scores()))
        format_generators.append(GenotypeGenerator(
            self._logger, self._gtc.get_gender(), self._gtc.get_genotypes()))
        return sorted(format_generators, cmp=self.format_cmp)

    def add_sample_format(self, bpm_records, vcf_record, sample_name):
        """
        Update the VCF record with the call
        data from a new sample

        Args:
            bpm_records (list(BPMRecord)): List of BPM records for this site
            vcf_record (vcf._Record): VCF record for this site
            sample_name (string): Name of the sample

        Returns:
            None

        Raises:
            Exception - No GTC file given for this call factory
        """
        if not self._gtc:
            raise Exception(
                "GTC file not available- cannot add sample information.")

        values = []

        for format_object in self._format_generators:
            values.append(format_object.generate_sample_format_info(
                bpm_records, vcf_record, sample_name))

        data = self._format_namedtuple(*values)
        call = _Call(vcf_record, sample_name, data)

        if vcf_record.samples is None:
            vcf_record.samples = []
        vcf_record._sample_indexes[call.sample] = len(vcf_record.samples)
        vcf_record.samples.append(call)
