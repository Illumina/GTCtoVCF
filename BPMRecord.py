from IlluminaBeadArrayFiles import RefStrand

COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "D": "D", "I": "I"}


def reverse(sequence):
    """
    Reverse a string

    Args:
        sequence (string): Sequence to reverse

    Returns:
        string: The reversed sequence
    """
    return sequence[::-1]


def complement(sequence):
    """
    Complement a nucleotide sequence. Note that complement of D and I are D and I, 
    respectively. This is intended to be called on the "SNP" portion of a source sequence.

    Args:
        sequence (string): The input sequence

    Returns:
        string: The complemented sequence
    """
    return "".join(COMPLEMENT_MAP[x] for x in sequence)


def reverse_complement(sequence):
    """
    Reverse complement a sequence

    Args:
        sequence (string): The input sequence

    Returns:
        string: The reverse-complement of the input sequence
    """
    return reverse(complement(sequence))


class BPMRecord(object):
    """
    Represents entry from a manifest.abs

    Attributes:
        name (string): Entry name (unique)
        address_a (string): Address of probe A
        chromosome (string): Chromsome name
        pos (int): Mapping position
        snp (string): SNP variation (e.g., [A/C])
        ref_strand (RefStrand): Reference strand of snp
        assay_type (int): 0 for for Inf II, 1 for Inf I
        indel_sequence (string): Sequence of indel, None for SNV
        index_num (int): Index in original manifest
        is_deletion (bool): Whether indel record represents deletion, None for SNV
    """

    def __init__(self, name, address_a, chromosome, pos, snp, ref_strand, assay_type, indel_sequence, index):
        self.name = name
        self.address_a = address_a
        self.chromosome = chromosome
        self.pos = int(pos)
        self.snp = snp
        self.ref_strand = ref_strand
        self.assay_type = assay_type
        self.indel_sequence = indel_sequence
        self.index_num = index
        self.is_deletion = None
        self._set_plus_strand_alleles()

    def is_indel(self):
        """
        Check whether a BPM record represents an indel

        Args:
            None

        Returns:
            bool: True if record is indel
        """
        return "D" in self.snp

    def get_plus_strand_indel_sequence(self):
        """
        Get the indel sequence on the plus strand

        Args:
            None

        Returns:
            string: Sequence of indel on plus strand. None for SNV
        """
        if self.indel_sequence is None:
            return None
        return self.indel_sequence if self.ref_strand == RefStrand.Plus else reverse_complement(self.indel_sequence)

    def _set_plus_strand_alleles(self):
        """
        Set the nucleotides alleles for the record on the plus strand.
        If record is indel, will set alleles in terms of D/I SNP convention
        
        Args:
            None
        
        Returns
            None
        
        Raises:
            Exception - Record does not contains reference strand information
        """
        nucleotides = [self.snp[1], self.snp[-2]]
        if self.ref_strand == RefStrand.Plus:
            self.plus_strand_alleles = nucleotides
        elif self.ref_strand == RefStrand.Minus:
            self.plus_strand_alleles = [
                COMPLEMENT_MAP[nucleotide] for nucleotide in nucleotides]
        else:
            raise Exception(
                "Manifest must contain reference strand information")

    def get_plus_strand_alleles(self):
        """
        Get the alleles for this record on the plus strand
        """
        return self.plus_strand_alleles
