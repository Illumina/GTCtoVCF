import os
from pysam import Fastafile


class ReferenceGenome(object):
    """
    Class to read sequence data from a reference genome

    Attributes:
        genome_fasta_file (string): Path to reference genome file
        fasta_file (pysam.Fastafile): Fastafile object for reference genome
    """

    def __init__(self, genome_fasta_file, logger):
        """
        Create new ReferenceGenome

        Args:
            genome_fasta_file (string): Path to whole genome FASTA file
            logger (logging.Logger): Logger for reporting warnings/errors

        Returns:
            ReferenceGenome
        """
        self._logger = logger
        self.genome_fasta_file = genome_fasta_file
        self._validate_reference_file(genome_fasta_file)

        try:
            self.fasta_file = Fastafile(self.genome_fasta_file)
        except:
            raise IOError("Could not read genome file: " +
                          self.genome_fasta_file)

    def _validate_reference_file(self, genome_fasta_file):
        """
        Check whether reference file is valid

        Args:
            genome_fasta_file (string): Path to genome FASTA file

        Returns
            bool: True if valid
        """
        if not os.path.isfile(self.genome_fasta_file + ".fai"):
            raise ValueError(
                "Supplied genome FASTA file does not have FASTA index")

    def get_reference_bases(self, chrom, start, end):
        """
        Get the reference bases from start to end

        Args:
            chrom (string): Chromsome to query
            start (int): Start position to query
            end (int): End position (not inclusive)

        Returns:
            string: The genome sequence

        Raises:
            ValueError - Invalid arguments
        """

        if start >= end:
            raise ValueError("Start/stop coordinates incorrect for: " +
                             str(chrom) + ":" + str(start) + "-" + str(end))

        if chrom not in self.fasta_file:
            raise ValueError(
                "FASTA reference is missing entry for chromosome " + str(chrom))

        return self.fasta_file.fetch(str(chrom), start, end)
