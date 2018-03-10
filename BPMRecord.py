from IlluminaBeadArrayFiles import RefStrand

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

    def __init__(self, name, address_a, probe_a, chromosome, pos, snp, ref_strand, assay_type, indel_source_sequence, source_strand, ilmn_strand, genome_reader, index):
        self.name = name
        self.address_a = address_a
        self.probe_a = probe_a
        self.chromosome = chromosome
        self.pos = int(pos)
        self.snp = snp
        self.ref_strand = ref_strand
        self.assay_type = assay_type
        self.indel_source_sequence = indel_source_sequence
        self._genome_reader = genome_reader
        self.index_num = index

        self.plus_strand_alleles = self._determine_plus_strand_alleles(snp, ref_strand)

        if self.indel_source_sequence:
            source_strand = source_strand[0].upper()
            ilmn_strand = ilmn_strand[0].upper()

            if source_strand == "U" or ilmn_strand == "U":
                raise ValueError("Unable to process indel with customer or ILMN strand value of \"U\"")

            if source_strand == "P" or source_strand == "M":
                assert ilmn_strand == "P" or ilmn_strand == "M"
            else:
                assert ilmn_strand == "T" or ilmn_strand == "B"

            self.is_source_on_design_strand = source_strand == ilmn_strand
            self.is_deletion = self._calculate_is_deletion()
        else:
            self.is_source_on_design_strand = None
            self.is_deletion = None

    def is_indel(self):
        """
        Check whether a BPM record represents an indel

        Args:
            None

        Returns:
            bool: True if record is indel
        """
        return "D" in self.snp

    def get_indel_source_sequences(self, ref_strand):
        return self.indel_source_sequence.get_split_sequence(self.is_source_on_design_strand != (self.ref_strand == ref_strand), True)

    def _calculate_is_deletion(self):
        if self.chromosome == "0" or self.pos == 0:
            return None

        start_index = self.pos - 1
        chromosome = "X" if self.chromosome == "XY" else self.chromosome

        # get indel sequence on the plus strand
        (five_prime, indel_sequence, three_prime) = self.get_indel_source_sequences(RefStrand.Plus)

        genomic_sequence = self._genome_reader.get_reference_bases(chromosome, start_index, start_index + len(indel_sequence))
        is_deletion = indel_sequence == genomic_sequence

        genomic_sequence = self._genome_reader.get_reference_bases(
            chromosome, start_index, start_index + len(indel_sequence))
        assert len(indel_sequence) == len(genomic_sequence)
        is_deletion = indel_sequence == genomic_sequence

        genomic_deletion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime), start_index)
        genomic_deletion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + len(indel_sequence), start_index + len(indel_sequence) + len(three_prime))

        genomic_insertion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime) + 1, start_index + 1)
        genomic_insertion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + 1, start_index + len(three_prime) + 1)

        deletion_context = max_suffix_match(genomic_deletion_five_prime, five_prime) + max_prefix_match(genomic_deletion_three_prime, three_prime)
        insertion_context = max_suffix_match(genomic_insertion_five_prime, five_prime) + max_prefix_match(genomic_insertion_three_prime, three_prime)

        if is_deletion and deletion_context > insertion_context:
            is_deletion = True
        elif insertion_context > deletion_context:
            is_deletion = False
        else:
            raise Exception("Unable to determine reference allele for indel")
        return is_deletion

    def _determine_plus_strand_alleles(self, snp, ref_strand):
        """
        Return the nucleotides alleles for the record on the plus strand.
        If record is indel, will return alleles in terms of D/I SNP convention

        Args:
            None

        Returns
            None

        Raises:
            Exception - Record does not contains reference strand information
        """
        nucleotides = [snp[1], snp[-2]]
        if ref_strand == RefStrand.Plus:
            return nucleotides
        elif ref_strand == RefStrand.Minus:
            return [
                COMPLEMENT_MAP[nucleotide] for nucleotide in nucleotides]
        else:
            raise Exception(
                "Manifest must contain reference strand information")

class IndelSourceSequence(object):
    """
    Represents the source sequence for an indel

    Attributes:
        five_prime (string) : Sequence 5' of indel (on the design strand)
        indel (string) : Indel sequence (on the design strand)
        three_prime (string) : Sequence 3' of indel (on the design strand)
    """
    def __init__(self, source_sequence):
        (self.five_prime, self.indel, self.three_prime) = self.split_source_sequence(source_sequence.upper())

    def get_split_sequence(self, generate_reverse_complement, left_shift):
        """
        Return the components of the indel source sequence

        Args:
            generate_reverse_complement (bool) : Return reverse complement of original source sequence
            left_shift (bool) : Left shift position of indel on requested strand

        Returns:
            (five_prime, indel, three_prime) = Sequences of three components of indel source sequence
        """
        if generate_reverse_complement:
            (five_prime, indel, three_prime) = (reverse_complement(self.three_prime), reverse_complement(self.indel), reverse_complement(self.five_prime))
        else:
            (five_prime, indel, three_prime) = (self.five_prime, self.indel, self.three_prime)

        if left_shift:
            while five_prime.endswith(indel):
                five_prime = five_prime[:-len(indel)]
                three_prime = indel + three_prime

        return (five_prime, indel, three_prime)

    @staticmethod
    def split_source_sequence(source_sequence):
        """
        Break source sequence into different piecdes

        Args:
            source_sequence (string): Source sequence string (e.g., ACGT[-/AGA]ATAT)

        Returns:
            (string, string, string) : Tuple with 5' sequence, indel sequence, 3' sequence
        """
        left_position = source_sequence.find("/")
        right_position = source_sequence.find("]")
        assert source_sequence[left_position - 1] == "-"
        return (source_sequence[:(left_position - 2)], source_sequence[(left_position + 1):right_position], source_sequence[(right_position+1):])

DEGENERACY_MAP = {}
DEGENERACY_MAP["A"] = "A"
DEGENERACY_MAP["C"] = "C"
DEGENERACY_MAP["G"] = "G"
DEGENERACY_MAP["T"] = "T"
DEGENERACY_MAP["R"] = "AG"
DEGENERACY_MAP["Y"] = "CT"
DEGENERACY_MAP["S"] = "GC"
DEGENERACY_MAP["W"] = "AT"
DEGENERACY_MAP["K"] = "GT"
DEGENERACY_MAP["M"] = "AC"
DEGENERACY_MAP["B"] = "CGT"
DEGENERACY_MAP["D"] = "AGT"
DEGENERACY_MAP["H"] = "ACT"
DEGENERACY_MAP["V"] = "ACG"
DEGENERACY_MAP["N"] = "ACGT"

def max_suffix_match(str1, str2):
    """
    Determine the maximum length of exact suffix
    match between str1 and str2

    str2 may contain degenerate IUPAC characters

    Args:
        str1 (string) : First string
        str2 (string) : Second string

    Returns:
        int : Length of maximum suffix match
    """
    result = 0
    for (char1, char2) in zip(str1[::-1], str2[::-1]):
        assert char1 in "ACGT"
        if char1 in DEGENERACY_MAP[char2]:
            result += 1
        else:
            break
    return result

def max_prefix_match(str1, str2):
    """
    Determine the maximum length of exact prefix
    match between str1 and str2

    str2 may contain degenerate IUPAC characters

    Args:
        str1 (string) : First string
        str2 (string) : Second string

    Returns:
        int : Length of maximum prefix match
    """
    result = 0
    for (char1, char2) in zip(str1, str2):
        assert char1 in "ACGT"
        if char1 in DEGENERACY_MAP[char2]:
            result += 1
        else:
            break
    return result
