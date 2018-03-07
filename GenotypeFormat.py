from itertools import combinations_with_replacement
from vcf.parser import _Format
from IlluminaBeadArrayFiles import code2genotype

CHANNEL_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def get_expected_ploidy(gender, chrom):
    """
    Determine expected ploidy of call based on sample's gender and chromosome. Unknown genders are processed as diploid.
        Args:
            gender (string): M,F or U
            chrom (string): chromosome, PAR values should be represented as XY
        Returns
            int: value of expected ploidy, currently set at 1 or 2
    """
    if (gender == "M" and chrom == "X") or chrom == "Y" or chrom == "MT":
        return 1

    return 2


def format_vcf_genotype(vcf_allele1_char, vcf_allele2_char, ploidy):
    """
    Create a VCF representation of the genotype based on alleles. Format appropriately if haploid.

    Args:
        vcf_allele1_char (string): 0,1,2 etc.
        vcf_allele2_char (string): 0,1,2, etc.
        vcf_record (vcf._Record): Record for the entry analyzed
        ploidy(int): Expected ploidy. 

    Returns
        string: String representation of genotype (e.g., "0/1"), adjusted for haploid calls if applicable
    """
    assert ploidy < 3

    # haploid calls
    if ploidy == 1:
        if vcf_allele1_char == vcf_allele2_char:
            return str(vcf_allele1_char)

    vcf_genotype = ""
    if vcf_allele2_char < vcf_allele1_char:
        vcf_genotype = str(vcf_allele2_char) + "/" + str(vcf_allele1_char)
    else:
        vcf_genotype = str(vcf_allele1_char) + "/" + str(vcf_allele2_char)
    return vcf_genotype


def convert_ab_genotype_to_nucleotide(ab_genotype, plus_strand_alleles):
    """
    Convert an integer genotype to nucleotide genotype on the plus strand

    Args:
        ab_genotype (int): 0 (no call), 1 (homozygous A), 2 (het), or 3 (homozygous B)
        plus_strand_alleles (tuple(string)): Tuple of length 2 with assay alleles (e.g., from SNP column of manifest)

    Returns
            tuple(string): Tuple of length 2 with nucleotide alleles on plus strand (e.g., ('A', 'C')). NC is ('-', '-')
    """
    if ab_genotype == 0:
        return ('-', '-')
    return tuple([plus_strand_alleles[0] if ab_allele == "A" else plus_strand_alleles[1] for ab_allele in code2genotype[ab_genotype]])


def convert_indel_genotype_to_vcf(nucleotide_genotypes, vcf_record, is_deletion, ploidy):
    """
    For indel, convert indel SNP genotype (e.g., I/D) into VCF genotype (e.g, 0/1)

    Args:
        nucleotide_genotypes (string,string): SNP genotype from manifest (e.g., ('D', 'I'))
        vcf_record (vcf._Record): Corresponding VCF record (define reference and alternate allele)
        is_deletion (bool): Whether the BPM record that produced the nucleotide genotypes is a reference deletion
        ploidy (int): Expected ploidy

    Returns:
        string: VCF genotype (e.g, "0/1")
    """
    if len(nucleotide_genotypes) == 0:
        return format_vcf_genotype(".", ".", ploidy)

    if len(nucleotide_genotypes) > 1:
        for nucleotide_genotype in nucleotide_genotypes:
            if nucleotide_genotype != nucleotide_genotypes[0]:
                return format_vcf_genotype(".", ".", ploidy)

    nucleotide_genotype = nucleotide_genotypes[0]
    if is_deletion:
        vcf_allele1_char = "0" if nucleotide_genotype[0] == "I" else "1"
        vcf_allele2_char = "0" if nucleotide_genotype[1] == "I" else "1"
    else:
        vcf_allele1_char = "1" if nucleotide_genotype[0] == "I" else "0"
        vcf_allele2_char = "1" if nucleotide_genotype[1] == "I" else "0"

    vcf_genotype = format_vcf_genotype(
        vcf_allele1_char, vcf_allele2_char, ploidy)

    return vcf_genotype


def convert_nucleotide_genotype_to_vcf(nucleotide_genotype, vcf_record, ploidy):
    """
    Convert a nucleotide genotype (on the plus strand) to a vcf genotype. For indels,
    use "convert_indel_genotype_to_vcf"

    Args:
        nucleotide_genotype (tuple(string)): Tuple of length 2. Each element is nucleotide allele
            No call should be represented as ("-", "-")
        ploidy(int): value of estimated ploidy.

    Returns
        string: VCF genotype (e.g, "0/1")
    """

    assert len(nucleotide_genotype) == 2
    if nucleotide_genotype[0] == "-" and nucleotide_genotype[1] == "-":
        return format_vcf_genotype(".", ".", ploidy)

    try:
        vcf_allele1_char = vcf_record.alleles.index(nucleotide_genotype[0])
        vcf_allele2_char = vcf_record.alleles.index(nucleotide_genotype[1])
    except:
        raise Exception(
            "Could not index alleles in VCF record " + vcf_record.ID)

    vcf_genotype = format_vcf_genotype(
        vcf_allele1_char, vcf_allele2_char, ploidy)

    return vcf_genotype


class RecordCombiner(object):
    """
    Class to take in a group of BPM records and output a combined genotype
    """

    def __init__(self, bpm_records, genotypes, logger):
        """
        Create new RecordCombiner

        Args:
            bpm_records (list(BPMRecord)): Group of BPM records for a single site (typically just one)
            genotypes (list(genotypes)): List of all genotypes in GTC file as integers

        Returns:
            RecordCombiner
        """
        self._bpm_records = bpm_records
        self._genotypes = genotypes
        self._logger = logger

    def _generate_possible_genotypes(self):
        """
        From the alleles in the BPM records, enumerate all possible genotypes
        at this site

        Args:
            None

        Returns:
            list(list(string)) - A list of lists of length 2. Each inner list represents a possible
                genotype at this site in terms of pair of nucleotide alleles on the plus strand
        """
        alleles = set()
        for record in self._bpm_records:
            alleles.update(record.get_plus_strand_alleles())
        return list(combinations_with_replacement(alleles, 2))

    def _record_inconsistent_with_genotype(self, record, genotype):
        """Check if a particular BPM record is inconsitent with a given
        genotype. Genotype should be tuple of alleles on plus strand.

        Args:
            record (BPMRecord)
            genotype (tuple(string)) : Tuple of length 2 where each element is nucleotide allele
                on plus strand
        Returns:
            bool
        """
        # record_genotype is an integer (0 - NC, 1 - AA, 2 - AB, 3 - BB)
        record_int_genotype = self._genotypes[record.index_num]
        if record_int_genotype == 0:
            return False

        plus_strand_alleles = record.get_plus_strand_alleles()
        record_plus_genotype = convert_ab_genotype_to_nucleotide(
            record_int_genotype, plus_strand_alleles)

        for allele in record_plus_genotype:
            # check for alleles that must be present in a consistent genotype
            consistent_alleles = []
            consistent_alleles.append(allele)
            if record.assay_type == 0:  # Inf II
                consistent_alleles.append(CHANNEL_MAP[allele])
            if not any([consistent_allele in genotype for consistent_allele in consistent_alleles]):
                return True

            # check for alleles that must be absent in a consistent genotype
            absent_alleles = []
            if record_int_genotype == 1 or record_int_genotype == 3:  # homozygous
                # for example, if assay is A/C Inf II and genotype is AA, real
                # genotype can not contain C or G
                absent_allele = plus_strand_alleles[0] if record_int_genotype == 3 else plus_strand_alleles[1]
                absent_alleles.append(absent_allele)
                if record.assay_type == 0:  # Inf II
                    absent_alleles.append(CHANNEL_MAP[absent_allele])
                if any([absent_allele in genotype for absent_allele in absent_alleles]):
                    return True
        return False

    def _filter_inconsistent_genotypes(self, possible_genotypes):
        """Filter the list of possible genotypes to remove
        those that are inconsitent with any BPM record in this group

        Args
            possible_genotypes (list(list(string)) - List of possible genotypes. Each possible genotype is a list of
                length 2 where each string is a nucleotide on the plus strand

        Returns
            list(list(string)) - List of genotypes consitent with assay data. Each remaining genotype is a list of
                length 2 where each string is a nucleotide on the plus strand
        """
        idx2inconsistent = [False] * len(possible_genotypes)
        for idx in xrange(len(possible_genotypes)):
            for record in self._bpm_records:
                if self._record_inconsistent_with_genotype(record, possible_genotypes[idx]):
                    idx2inconsistent[idx] = True
                    break
        return [genotype for (genotype, is_inconsistent) in zip(possible_genotypes, idx2inconsistent) if not is_inconsistent]

    def combine_genotypes(self):
        """
        Generate the combined genotype from all assays at this site

        Args:
            None

        Returns:
            (string, string): The combined genotype (on the plus strand) at this site (e.g., ("A", "C") ) No call is ("-", "-")
        """
        possible_genotypes = self._generate_possible_genotypes()

        # filter any genotypes which may be ambiguous due to the presence of
        # InfII assays
        allowable_genotypes = self._filter_inconsistent_genotypes(
            possible_genotypes)

        # if only one consistent genotype, then that is the genotype. Otherwise, the genotype
        # is ambiguous (more than 1) or inconsistent (less than 1) and return a
        # no-call
        if len(allowable_genotypes) == 1:
            return allowable_genotypes[0]
        return ("-", "-")

    def combine_names(self):
        """
        Generate the combined name for thi sgroup of records

        Args:
            None

        Returns:
            string: The combined names
        """
        record_names = []
        for record in self._bpm_records:
            record_names.append(record.name)
        return ",".join(sorted(record_names))


class GenotypeFormat(object):
    """
    Generate GT format information for VCF
    """
    def __init__(self, logger, gender, genotypes):
        self._gender = gender
        self._genotypes = genotypes
        self._logger = logger

    @staticmethod
    def get_id():
        return "GT"

    @staticmethod
    def get_description():
        return "Genotype"

    @staticmethod
    def get_format_obj():
        return _Format(GenotypeFormat.get_id(), 1, "String", GenotypeFormat.get_description())

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """ 
        Get the sample genotype

        Args:
            bpm_records (list(BPMRecord)): List of BPM records
            vcf_record (vcf._Record): Corresponding VCF record
            sample_name (string): The sample name

        Returns:
            string: GT sample format string (e.g., "0/1")
        """
        # use chrom from bpm_records where PAR chromsomes have not been
        # converted
        ploidy = get_expected_ploidy(self._gender, bpm_records[0].chromosome)
        if any(record.is_indel() for record in bpm_records):
            assert all(record.is_indel() for record in bpm_records)
            nucleotide_genotypes = []
            for record in bpm_records:
                int_genotype = self._genotypes[record.index_num]
                if int_genotype != 0:
                    nucleotide_genotypes.append(convert_ab_genotype_to_nucleotide(
                        int_genotype, bpm_records[0].get_plus_strand_alleles()))
            vcf_genotype = convert_indel_genotype_to_vcf(
                nucleotide_genotypes, vcf_record, bpm_records[0].is_deletion, ploidy)
        else:
            if len(bpm_records) > 1:
                combiner = RecordCombiner(
                    bpm_records, self._genotypes, self._logger)
                nucleotide_genotype = combiner.combine_genotypes()
                vcf_genotype = convert_nucleotide_genotype_to_vcf(
                    nucleotide_genotype, vcf_record, ploidy)
                vcf_record.ID = combiner.combine_names()
            else:
                sample_genotype = self._genotypes[bpm_records[0].index_num]
                if sample_genotype == 0:
                    nucleotide_genotype = ('-', '-')
                else:
                    nucleotide_genotype = convert_ab_genotype_to_nucleotide(
                        sample_genotype, bpm_records[0].get_plus_strand_alleles())
                vcf_genotype = convert_nucleotide_genotype_to_vcf(
                    nucleotide_genotype, vcf_record, ploidy)

        return vcf_genotype
