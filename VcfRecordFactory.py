from IlluminaBeadArrayFiles import RefStrand
from vcf.model import _Record, _Substitution

def check_reference_allele(reference_base, bpm_record_group):
    """
    Check whether the given reference base (on the plus strand) is queried
    by any record in a group of BPMRecords

    Args:
        reference_base (string): The reference base
        bpm_record_group (iter(BPMRecord)): Iterable of BPMRecords

    Returns:
        bool: True if record group contains reference allele
    """
    if any([reference_base in record.get_plus_strand_alleles() for record in bpm_record_group]):
        return True
    return False


class VcfRecordFactory(object):
    """Class to create new VCF records"""

    def __init__(self, format_factory, genome_reader, expand_identifiers, auxiliary_records, logger):
        """
        Create a new VcfRecordFactory

        Args:
            format_factory (FormatFactory): Provides format information need to create VCF records
            genome_reader (GenomeReader): Query genomic sequence regions
            expand_identifiers (bool): If true, VCF record id will be comma-delimited list of all BPM identifiers. If false, will use shortest identifier.
            auxiliary_records (dict(string, vcf._Record)): Maps from BPM identifer to VCF record with alternate annotation for that BPM entry
            logger (logging.Logger): Logger to report warnings/errors

        Returns:
            VcfRecordFactory
        """
        self._format_factory = format_factory
        self._genome_reader = genome_reader
        self._expand_identifiers = expand_identifiers
        self._auxiliary_records = auxiliary_records
        self._logger = logger

    def create_vcf_record(self, bpm_record_group):
        """
        Generate an VCF record from a group of BPM records.
        Note that this may alter the BPM records for any
        auxiliary loci

        Args:
            bpm_record_group (list(BPMRecord)): List of BPM records to use to create single VCF record

        Returns:
            vcf._Record: New VCF record
        """
        bpm_record = bpm_record_group[0]

        qual = ""
        filt = "PASS"
        info = []
        sample_indexes = {}
        start_index = bpm_record.pos - 1
        chrom = bpm_record.chromosome

        if len(bpm_record_group) == 1:
            identifier = bpm_record.name
        elif self._expand_identifiers:
            identifier = ",".join(
                [bpm_record.name for bpm_record in bpm_record_group])
        else:
            min_identifier_length = min(
                [len(bpm_record.name) for bpm_record in bpm_record_group])
            identifier = [bpm_record.name for bpm_record in bpm_record_group if len(
                bpm_record.name) == min_identifier_length][0]

        # convert PAR annotations
        if chrom == "XX" or chrom == "XY":
            chrom = "X"

        if self._auxiliary_records is not None and any([bpm_record.name in self._auxiliary_records for bpm_record in bpm_record_group]):
            if len(bpm_record_group) > 1:
                raise Exception(
                    "Unable to use auxiliary definition of locus " + identifier + " with multiple assays")
            elif bpm_record.is_indel():
                raise Exception(
                    "Unable to use auxiliary definition of locus " + identifier + " with indel assay")
            else:
                auxiliary_record = self._auxiliary_records[bpm_record.name]
                if len(auxiliary_record.alleles) != 2:
                    raise Exception(
                        "Auxiliary locus definition for " + auxiliary_record.ID + " is not bi-allelic")
                for allele in auxiliary_record.alleles:
                    if len(str(allele)) <= 1:
                        raise Exception("Auxiliary locus definition for " +
                                        auxiliary_record.ID + " is not a multi-nucleotide variant")
                else:
                    # update BPM record with plus-strand alleles
                    # these must stay in the same order as so plus allele shouldn't be reported (necessarily) as reference, alternate
                    # need to use ref strand info to figure out which end of MNV allele to use for comparisonn
                    # then create vcf record with alleles from auxiliary record
                    old_plus_strand_alleles = bpm_record.get_plus_strand_alleles()
                    new_plus_strand_alleles = []
                    if bpm_record.ref_strand == RefStrand.Plus:
                        for old_allele in old_plus_strand_alleles:
                            for new_allele in auxiliary_record.alleles:
                                if str(new_allele)[0] == old_allele[0]:
                                    new_plus_strand_alleles.append(
                                        str(new_allele))
                                    break
                    else:
                        assert bpm_record.ref_strand == RefStrand.Minus
                        for old_allele in old_plus_strand_alleles:
                            for new_allele in auxiliary_record.alleles:
                                if str(new_allele)[-1] == old_allele[-1]:
                                    new_plus_strand_alleles.append(
                                        str(new_allele))
                                    break
                    assert len(new_plus_strand_alleles) == 2
                    assert new_plus_strand_alleles[0] != new_plus_strand_alleles[1]
                    bpm_record.plus_strand_alleles = new_plus_strand_alleles

                    return _Record(auxiliary_record.CHROM, auxiliary_record.POS, identifier, auxiliary_record.REF, auxiliary_record.ALT, qual, filt, info,
                                   self._format_factory.get_format_id_string(), sample_indexes)
        else:
            if bpm_record.is_indel():
                # note that indel will have been shifted five prime, such that five_prime sequence can not
                # end with indel_sequence
                (five_prime, indel_sequence, three_prime) = bpm_record.indel_source_sequence.get_plus_strand_sequence(bpm_record.ref_strand)
                for record in bpm_record_group:
                    assert record.is_indel() and record.indel_source_sequence.get_plus_strand_sequence(record.ref_strand)[1] == indel_sequence
                
                genomic_sequence = self._genome_reader.get_reference_bases(
                    chrom, start_index, start_index + len(indel_sequence))
                assert len(indel_sequence) == len(genomic_sequence)
                is_deletion = indel_sequence == genomic_sequence

                def max_suffix_match(str1, str2):
                    result = 0
                    for (x,y) in zip(str1[::-1], str2[::-1]):
                        if x == y:
                            result += 1
                        else:
                            break
                    return result
                
                def max_prefix_match(str1, str2):
                    result = 0
                    for (x,y) in zip(str1, str2):
                        if x == y:
                            result += 1
                        else:
                            break
                    return result

                genomic_deletion_five_prime = self._genome_reader.get_reference_bases(
                    chrom, start_index - len(five_prime), start_index)
                genomic_deletion_three_prime = self._genome_reader.get_reference_bases(
                    chrom, start_index + len(indel_sequence), start_index + len(indel_sequence) + len(three_prime))

                genomic_insertion_five_prime = self._genome_reader.get_reference_bases(
                    chrom, start_index - len(five_prime) + 1, start_index + 1)
                genomic_insertion_three_prime = self._genome_reader.get_reference_bases(
                    chrom, start_index + 1, start_index + len(three_prime) + 1)

                deletion_context = max_suffix_match(genomic_deletion_five_prime, five_prime) + max_prefix_match(genomic_deletion_three_prime, three_prime)
                insertion_context = max_suffix_match(genomic_insertion_five_prime, five_prime) + max_prefix_match(genomic_insertion_three_prime, three_prime)

                if is_deletion and deletion_context > insertion_context:
                    is_deletion = True
                elif insertion_context > deletion_context:
                    is_deletion = False
                else:
                    raise Exception("Unable to determine reference allele for indel")                

                for record in bpm_record_group:
                    record.is_deletion = is_deletion

                if is_deletion:
                    reference_base = self._genome_reader.get_reference_bases(
                        chrom, start_index - 1, start_index)
                    reference_allele = reference_base + indel_sequence
                    alternate_allele = reference_base
                    return _Record(chrom, start_index, identifier, reference_allele,
                                   [_Substitution(alternate_allele)
                                    ], qual, filt, info,
                                   self._format_factory.get_format_id_string(), sample_indexes)
                else:
                    reference_base = self._genome_reader.get_reference_bases(
                        chrom, start_index, start_index + 1)
                    reference_allele = reference_base
                    alternate_allele = reference_base + indel_sequence
                    return _Record(chrom, start_index + 1, identifier, reference_allele,
                                   [_Substitution(alternate_allele)
                                    ], qual, filt, info,
                                   self._format_factory.get_format_id_string(), sample_indexes)
            else:
                reference_base = self._genome_reader.get_reference_bases(
                    chrom, start_index, start_index + 1)
                if not check_reference_allele(reference_base, bpm_record_group):
                    self._logger.warn(
                        "Reference allele is not queried for locus: " + identifier)

                alts = []
                for record in bpm_record_group:
                    for nucleotide in record.get_plus_strand_alleles():
                        if nucleotide != reference_base:
                            substitution = _Substitution(nucleotide)
                            if substitution not in alts:
                                alts.append(_Substitution(nucleotide))

                return _Record(chrom, bpm_record.pos, identifier, reference_base, alts, qual, filt, info,
                               self._format_factory.get_format_id_string(), sample_indexes)
