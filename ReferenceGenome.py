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
        self._validate_reference_file()

        try:
            self._fasta_file = Fastafile(self.genome_fasta_file)
        except:
            raise IOError("Could not read genome file: " +
                          self.genome_fasta_file)

    def _validate_reference_file(self):
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

    def get_contig_lengths(self):
        """
        Get the names and lengths of all contigs in a references

        Args:
            None

        Returns:
            list(tuple(string,int)): Returns a list representing the name and lengths of all contigs in reference
        """
        return zip(self._fasta_file.references, self._fasta_file.lengths)

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

        if chrom not in self._fasta_file:
            raise ValueError(
                "FASTA reference is missing entry for chromosome " + str(chrom))

        return self._fasta_file.fetch(str(chrom), start, end)

class CachedReferenceGenome(object):
    """
    Class to provide sequence data from a reference genome. All sequence
    data will be cached into memory.

    Attributes:
        genome_fasta_file (string): Path to reference genome file
    """
    def __init__(self, reference_genome, logger): 
        self._logger = logger
        self._logger.info("Caching reference data")
        self._cache = CachedReferenceGenome.generate_genome_cache(reference_genome)
        self._logger.info("Finished caching reference data")
        self.genome_fasta_file = reference_genome.genome_fasta_file

    @staticmethod
    def generate_genome_cache(reference_genome):
        """
        Get the reference bases from start to end

        Args:
            reference_genome (ReferenceGenome): A reference genome that provides the information to be cached

        Returns:
            dict(string,string): A dictionary that maps from contig name to bases for that contig
        """        
        result = {}
        for (contig, contig_length) in reference_genome.get_contig_lengths():
            result[contig] = reference_genome.get_reference_bases(contig, 0, contig_length)
        return result 

    def get_contig_lengths(self):
        """
        Get the names and lengths of all contigs in a references

        Args:
            None

        Returns:
            list(tuple(string,int)): Returns a list representing the name and lengths of all contigs in reference
        """
        return [(contig, len(self._cache[contig])) for contig in self._cache]

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

        if chrom not in self._cache:
            raise ValueError(
                "Reference is missing entry for chromosome " + str(chrom))

        return self._cache[str(chrom)][start:end]
