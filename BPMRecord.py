from IlluminaBeadArrayFiles import RefStrand

#COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "D": "D", "I": "I"}

COMPLEMENT_MAP = dict(zip("ABCDGHKMRTVYNID", "TVGHCDMKYABRNID"))

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
        probe_a (string): Sequence of probe A
        chromosome (string): Chromsome name
        pos (int): Mapping position
        snp (string): SNP variation (e.g., [A/C])
        ref_strand (RefStrand): Reference strand of snp
        assay_type (int): 0 for for Inf II, 1 for Inf I
        indel_source_sequence (IndelSourceSequence): Sequence of indel (on design strand), None for SNV
        index_num (int): Index in original manifest
        is_deletion (bool): Whether indel record represents deletion, None for SNV
    """

    def __init__(self, name, address_a, probe_a, chromosome, pos, snp, ref_strand, assay_type, indel_source_sequence, index):
        self.name = name
        self.address_a = address_a
        self.probe_a = probe_a
        self.chromosome = chromosome
        self.pos = int(pos)
        self.snp = snp
        self.ref_strand = ref_strand
        self.assay_type = assay_type
        self.indel_source_sequence = indel_source_sequence
        self.index_num = index
        self.is_deletion = None
        self._set_plus_strand_alleles()

    def get_plus_strand_probe(self, shorten = False):
        result = self.probe_a[:-1] if shorten and self.assay_type == 1 else self.probe_a
        return result if self.ref_strand == RefStrand.Plus else reverse_complement(result)

    def is_indel(self):
        """
        Check whether a BPM record represents an indel

        Args:
            None

        Returns:
            bool: True if record is indel
        """
        return "D" in self.snp

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


class IndelSourceSequence(object):
    """
    Represents the source sequence for an indel

    Attributes:
        five_prime (string) : Ssequence 5' of indel (on the design strand)
        indel (string) : Indel sequence (on the design strand)
        three_prime (string) : Sequence 3' of indel (on the design strand)
    """
    def __init__(self, source_sequence, source_strand, ilmn_strand):
        assert "-" in source_sequence
        source_sequence = source_sequence.upper()

        source_strand = source_strand[0].upper()
        ilmn_strand = ilmn_strand[0].upper()

        if source_strand == "U" or ilmn_strand == "U":
            raise ValueError("Unable to process indel with customer or ILMN strand value of \"U\"")

        (five_prime, indel, three_prime) = split_source_sequence(source_sequence)

        if source_strand == "P" or source_strand == "M":
            assert ilmn_strand == "P" or ilmn_strand == "M"
        else:
            assert ilmn_strand == "T" or ilmn_strand == "B"

        if source_strand != ilmn_strand:
            (self.five_prime, self.indel, self.three_prime) = (reverse_complement(three_prime), reverse_complement(indel), reverse_complement(five_prime))
        else:
            (self.five_prime, self.indel, self.three_prime) = (five_prime, indel, three_prime)

    def get_plus_strand_sequence(self, ref_strand):
        """
        Position will be left shifted (five_prime will not end with indel)
        """
        if ref_strand == RefStrand.Plus:
            (five_prime, indel, three_prime) = (self.five_prime, self.indel, self.three_prime)
        elif ref_strand == RefStrand.Minus:
            (five_prime, indel, three_prime) = (reverse_complement(self.three_prime), reverse_complement(self.indel), reverse_complement(self.five_prime))
        else:
            raise RuntimeError("Unknown ref_strand value " + ref_strand)
        
        while five_prime.endswith(indel):
            five_prime = five_prime[:-len(indel)]
            three_prime = indel + three_prime
        
        return (five_prime, indel, three_prime)

def split_source_sequence(source_sequence):
    """
    Break source sequence into different piecdes

    Args:
        source_sequence (string): Source sequence string

    Returns:
        (string, string, string) : Tuple with 5' sequence, indel sequence, 3' sequence
    """
    left_position = source_sequence.find("/")
    right_position = source_sequence.find("]")
    assert source_sequence[left_position - 1] == "-"
    return (source_sequence[:(left_position - 2)], source_sequence[(left_position + 1):right_position],source_sequence[(right_position+1):])