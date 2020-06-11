import os
from tempfile import mkstemp
from collections import OrderedDict
from vcf.parser import Reader, _Contig

class ReaderTemplateFactory(object):
    """Class to create new reader templates"""

    def __init__(self, genome_reader, format_factory, vcf_version, software_id, chrom_order, logger):
        """Create new ReaderTemplateFactory

        Args
            genome_reader (GenomeReader): Genome reader
            format_factory (FormatFactory): FormatFactory to provide information about format_string
            vcf_version (string): version of VCF file
            software_id (string): identifier for software creating the VCF
            chrom_order dict(string,int) : dictionary from contig to order
            chrom_sort_function(string): function definition used to sort chromosomes
        Returns
            ReaderTemplateFactory
        """
        self._genome_reader = genome_reader
        self._format_factory = format_factory
        self._vcf_version = vcf_version
        self._software_id = software_id
        self._chrom_order = chrom_order
        self._contigs = {}
        self._logger = logger

    def add_format_data(self, reader_template):
        """
        Query the CallFactory object to get FORMAT info

        Args:
            reader_template (vcf.Reader): Reader template to be populated with FORMAT info

        Returns:
            None
        """

        format_objects = self._format_factory.get_format_objects()
        if not format_objects:
            self._logger.warn("No format info available")
            return

        for format_object in format_objects:
            reader_template.formats[format_object.id] = format_object

    def create_reader_template(self, sample_ids):
        """
        Create a new reader template (vcf.Reader)

        Args:
            sample_ids (list(string)) : List of sample IDs for header of template

        Returns:
            Reader
        """

        file_descriptor, temp_path = mkstemp(suffix=".vcf")
        with os.fdopen(file_descriptor, "w") as vcf_handle:
            vcf_handle.write("##fileformat=VCFv" + self._vcf_version + "\n")
            vcf_handle.write("##source=" + self._software_id + "\n")
            vcf_handle.write(
                "##reference=file://" + os.path.abspath(self._genome_reader.genome_fasta_file) + "\n")
            vcf_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        with open(temp_path) as vcf_handle:
            reader_template = Reader(vcf_handle)

        reader_template.contigs = self.create_header_contigs()
        if sample_ids:
            reader_template._column_headers.append("FORMAT")
            self.add_format_data(reader_template)
            reader_template.samples = sample_ids
            reader_template._sample_indexes = dict(
                [(x, i) for (i, x) in enumerate(reader_template.samples)])

        os.remove(temp_path)

        return reader_template

    def create_header_contigs(self):
        """
        Create contig information for the header

        Args:
            None

        Returns
            dict(string,_Contig): Contig information
        """
        contigs = {}
        for contig_name, length in self._genome_reader.get_contig_lengths():
            if contig_name.isdigit():
                contigs[int(contig_name)] = _Contig(contig_name, length)
            else:
                contigs[contig_name] = _Contig(contig_name, length)
        contigs = OrderedDict(
            sorted(contigs.items(), key=lambda x: self._chrom_order[str(x[0])]))
        return contigs
