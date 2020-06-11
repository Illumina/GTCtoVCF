from .IlluminaBeadArrayFiles import RefStrand

COMPLEMENT_MAP = dict(zip("ABCDGHKMRTVYNID", "TVGHCDMKYABRNID"))
REQUIRED_INDEL_CONTEXT_LENGTH = 3


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


def determine_left_shift(five_prime, indel, three_prime):
    """
    Adjust 5' and 3' context of indel such that
    indel is fully shifted to 5'

    Args:
        five_prime (string) : Five prime sequence
        indel (string) : Sequence of indel
        three_prime (string) : Three prime sequence

    Returns:
        (string, string) : New sequence of 5' and 3' sequences
    """
    while five_prime.endswith(indel):
        five_prime = five_prime[:-len(indel)]
        three_prime = indel + three_prime
    # may have not fully shifted homopolymer
    while len(indel) * five_prime[-1] == indel:
        three_prime = five_prime[-1] + three_prime
        five_prime = five_prime[:-1]
    return five_prime, three_prime


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
        assay_type (int): 0 for Inf II, 1 for Inf I
        indel_source_sequence (IndelSourceSequence): Sequence of indel (on design strand), None for SNV
        index_num (int): Index in original manifest
        is_deletion (bool): Whether indel record represents deletion, None for SNV
    """

    def __init__(self, name, address_a, probe_a, chromosome, pos, snp, ref_strand, assay_type, indel_source_sequence,
                 source_strand, ilmn_strand, genome_reader, index, logger):
        """
        Create a new BPM record

        Args:
            name (string) : Name field from manifest
            address_a (string) : AddressA_ID field from manifest
            probe_a (string) : AlleleA_ProbeSeq field from manifest
            chromosome (string) : Chr field from manifest
            pos (string, int) : MapInfo field from manifest
            ref_strand (RefStrand) : RefStrand from manifest
            assay_type (int) : 0 for Inf II, 1 for Inf I
            indel_source_sequence (IndelSourceSequence) : Source sequence for indel, may be None for SNV
            source_strand (string) : SourceStrand field from manifest
            ilmn_strand (strinng) : IlmnStrand field from manifest
            genome_reader (ReferenceGenome,CachedReferenceGenome) : Allows query of genomic sequence, may be None for SNV
            index (int) : Index of entry within manifest/GTC files
            logger (Logger) : A logger
        """
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
        self._logger = logger

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

        Returns:
            bool: True if record is indel
        """
        return "D" in self.snp

    def get_indel_source_sequences(self, ref_strand):
        return self.indel_source_sequence.get_split_sequence(
            self.is_source_on_design_strand != (self.ref_strand == ref_strand), True)

    def _calculate_is_deletion(self):
        if self.chromosome == "0" or self.pos == 0:
            return None

        start_index = self.pos - 1
        chromosome = "X" if self.chromosome == "XY" else self.chromosome

        # get indel sequence on the plus strand
        (five_prime, indel_sequence, three_prime) = self.get_indel_source_sequences(RefStrand.Plus)

        genomic_sequence = self._genome_reader.get_reference_bases(chromosome, start_index,
                                                                   start_index + len(indel_sequence))
        indel_sequence_match = indel_sequence == genomic_sequence

        genomic_deletion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime), start_index)
        genomic_deletion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + len(indel_sequence), start_index + len(indel_sequence) + len(three_prime))
        (genomic_deletion_five_prime, genomic_deletion_three_prime) = determine_left_shift(genomic_deletion_five_prime,
                                                                                           indel_sequence,
                                                                                           genomic_deletion_three_prime)

        genomic_insertion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime) + 1, start_index + 1)
        genomic_insertion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + 1, start_index + len(three_prime) + 1)
        (genomic_insertion_five_prime, genomic_insertion_three_prime) = determine_left_shift(
            genomic_insertion_five_prime, indel_sequence, genomic_insertion_three_prime)

        deletion_context_match_lengths = (max_suffix_match(genomic_deletion_five_prime, five_prime),
                                          max_prefix_match(genomic_deletion_three_prime, three_prime))
        max_deletion_context = min(len(genomic_deletion_five_prime), len(five_prime)) + min(
            len(genomic_deletion_three_prime), len(three_prime)) + len(indel_sequence)
        deletion_context_score = (sum(deletion_context_match_lengths) + len(
            indel_sequence) if indel_sequence_match else 0) / float(max_deletion_context)

        insertion_context_match_lengths = (max_suffix_match(genomic_insertion_five_prime, five_prime),
                                           max_prefix_match(genomic_insertion_three_prime, three_prime))
        max_insertion_context = min(len(genomic_insertion_five_prime), len(five_prime)) + min(
            len(genomic_insertion_three_prime), len(three_prime))
        insertion_context_score = sum(insertion_context_match_lengths) / float(max_insertion_context)

        is_deletion = indel_sequence_match and deletion_context_score > insertion_context_score and min(
            deletion_context_match_lengths) >= 1

        is_insertion = insertion_context_score > deletion_context_score and min(insertion_context_match_lengths) >= 1

        if is_deletion == is_insertion:
            raise Exception("Unable to determine reference allele for indel")

        if is_deletion:
            if deletion_context_score < 1.0:
                self._logger.warn("Incomplete match of source sequence to genome for indel " + self.name)

        if is_insertion:
            if insertion_context_score < 1.0:
                self._logger.warn("Incomplete match of source sequence to genome for indel " + self.name)

        return is_deletion

    @staticmethod
    def _determine_plus_strand_alleles(snp, ref_strand):
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
            (five_prime, indel, three_prime) = (
                reverse_complement(self.three_prime), reverse_complement(self.indel),
                reverse_complement(self.five_prime))
        else:
            (five_prime, indel, three_prime) = (self.five_prime, self.indel, self.three_prime)

        if left_shift:
            (five_prime, three_prime) = determine_left_shift(five_prime, indel, three_prime)
        return five_prime, indel, three_prime

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
        return (source_sequence[:(left_position - 2)], source_sequence[(left_position + 1):right_position],
                source_sequence[(right_position + 1):])


DEGENERACY_MAP = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
                  "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


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
