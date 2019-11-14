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
    if any([reference_base in record.plus_strand_alleles for record in bpm_record_group]):
        return True
    return False

class VcfRecordFactory(object):
    """
    Class to create new VCF records
    """

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
        Generate a VCF record from a group of BPM records.
        Note that this may alter the BPM records for any
        auxiliary loci

        Args:
            bpm_record_group (list(BPMRecord)): List of BPM records to use to create single VCF record

        Returns:
            vcf._Record: New VCF record
        """
        bpm_record = bpm_record_group[0]

        if self._auxiliary_records and any([bpm_record.name in self._auxiliary_records for bpm_record in bpm_record_group]):
            if len(bpm_record_group) > 1:
                raise Exception(
                    "Unable to use auxiliary definition of locus " + self._get_identifier(bpm_record_group) + " with multiple assays")
            elif bpm_record.is_indel():
                raise Exception(
                    "Unable to use auxiliary definition of locus " + self._get_identifier(bpm_record_group) + " with indel assay")
            else:
                auxiliary_record = self._auxiliary_records[bpm_record.name]
                if len(auxiliary_record.alleles) != 2:
                    raise Exception(
                        "Auxiliary locus definition for " + auxiliary_record.ID + " is not bi-allelic")
                if any([len(str(allele)) <= 1 for allele in auxiliary_record.alleles]):
                    raise Exception("Auxiliary locus definition for " +
                                        auxiliary_record.ID + " is not a multi-nucleotide variant")
                return self._get_record_for_auxiliary(bpm_record_group, auxiliary_record)
        else:
            if bpm_record.is_indel():
                return self._get_record_for_indel(bpm_record_group)
            return self._get_record_for_snv(bpm_record_group)

    def _get_identifier(self, bpm_record_group):
        """
        Generate an single identifier for a group of BPM record groups

        Args:
            bpm_record_group : The group of BPM record

        Returns:
            string : the identifier
        """
        if len(bpm_record_group) == 1:
            return bpm_record_group[0].name
        elif self._expand_identifiers:
            return ";".join(
                [bpm_record.name for bpm_record in bpm_record_group])
        min_identifier_length = min(
            [len(bpm_record.name) for bpm_record in bpm_record_group])
        return [bpm_record.name for bpm_record in bpm_record_group if len(
            bpm_record.name) == min_identifier_length][0]

    @staticmethod
    def _get_record_defaults():
        """
        Generate some default elements for the VCF record

        Args:
            None

        Returns:
            (string, string, list, dict) : ("", "PASS", [], {})
        """
        return ("", "PASS", [], {})

    def _get_record_for_auxiliary(self, bpm_record_group, auxiliary_record):
        """
        Create a new VCF record for an auxiliar locus
        THIS WILL ALSO ALTER THE BPM RECORD

        Args:
            bpm_record_group (list(BPMRecord)) : BPM records for the group (must be length 1)
            auxiliary_record (vcf._Record) : Auxiliar locus definition

        Returns:
            vcf._Record : The new VCF record definition
        """
        # update BPM record with plus-strand alleles
        # these must stay in the same order as so plus allele shouldn't be reported (necessarily) as reference, alternate
        # need to use ref strand info to figure out which end of MNV allele to use for comparisonn
        # then create vcf record with alleles from auxiliary record
        identifier = self._get_identifier(bpm_record_group)

        assert len(bpm_record_group) == 1
        bpm_record = bpm_record_group[0]
        (qual, filt, info, sample_indexes) = VcfRecordFactory._get_record_defaults()
        old_plus_strand_alleles = bpm_record.plus_strand_alleles
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

    def _get_record_for_snv(self, bpm_record_group):
        """
        Create a new VCF record for an SNV

        Args:
            bpm_record_group (list(BPMRecord)) : BPM records for the group (must all be SNV records)

        Returns:
            vcf._Record : The new VCF record definition
        """
        assert not bpm_record_group[0].is_indel()
        identifier = self._get_identifier(bpm_record_group)
        bpm_record = bpm_record_group[0]
        (qual, filt, info, sample_indexes) = VcfRecordFactory._get_record_defaults()

        start_index = bpm_record.pos - 1
        chrom = bpm_record.chromosome
        if chrom == "XX" or chrom == "XY":
            chrom = "X"

        reference_base = self._genome_reader.get_reference_bases(
            chrom, start_index, start_index + 1)
        if not check_reference_allele(reference_base, bpm_record_group):
            self._logger.warn(
                "Reference allele is not queried for locus: " + identifier)

        alts = []
        for record in bpm_record_group:
            for nucleotide in record.plus_strand_alleles:
                if nucleotide != reference_base:
                    substitution = _Substitution(nucleotide)
                    if substitution not in alts:
                        alts.append(_Substitution(nucleotide))

        return _Record(chrom, bpm_record.pos, identifier, reference_base, alts, qual, filt, info,
                       self._format_factory.get_format_id_string(), sample_indexes)

    def _get_record_for_indel(self, bpm_record_group):
        """
        Create a new VCF record for an indel

        Args:
            bpm_record_group (list(BPMRecord)) : BPM records for the group (must all be indels)

        Returns:
            vcf._Record : The new VCF record definition
        """
        (qual, filt, info, sample_indexes) = VcfRecordFactory._get_record_defaults()

        identifier = self._get_identifier(bpm_record_group)
        for record in bpm_record_group:
            assert record.is_deletion == bpm_record_group[0].is_deletion

        bpm_record = bpm_record_group[0]

        (_, indel_sequence, _) = bpm_record.get_indel_source_sequences(RefStrand.Plus)
        start_index = bpm_record.pos - 1
        chrom = bpm_record.chromosome
        if chrom == "XX" or chrom == "XY":
            chrom = "X"

        if bpm_record.is_deletion:
            reference_base = self._genome_reader.get_reference_bases(
                chrom, start_index - 1, start_index)
            reference_allele = reference_base + indel_sequence
            alternate_allele = reference_base
            return _Record(chrom, start_index, identifier, reference_allele,
                           [_Substitution(alternate_allele)
                            ], qual, filt, info,
                           self._format_factory.get_format_id_string(), sample_indexes)

        reference_base = self._genome_reader.get_reference_bases(
            chrom, start_index, start_index + 1)
        reference_allele = reference_base
        alternate_allele = reference_base + indel_sequence
        return _Record(chrom, start_index + 1, identifier, reference_allele,
                       [_Substitution(alternate_allele)
                        ], qual, filt, info,
                       self._format_factory.get_format_id_string(), sample_indexes)
