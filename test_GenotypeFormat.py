import unittest
from GenotypeFormat import convert_indel_genotype_to_vcf

class DummyRecord:
    # convert_indel_genotype_to_vcf does not use the vcf_record, but keep a placeholder
    pass

class TestConvertIndelGenotypeToVcf(unittest.TestCase):
    def test_mixed_DI_ID_returns_nocall_diploid_deletion_true(self):
        nucleotide_genotypes = [('D', 'I'), ('I', 'D')]
        result = convert_indel_genotype_to_vcf(nucleotide_genotypes, DummyRecord(), True, 2)
        self.assertEqual(result, "./.")

    def test_mixed_DI_ID_returns_nocall_diploid_deletion_false(self):
        nucleotide_genotypes = [('D', 'I'), ('I', 'D')]
        result = convert_indel_genotype_to_vcf(nucleotide_genotypes, DummyRecord(), False, 2)
        self.assertEqual(result, "./.")

    def test_mixed_DI_ID_returns_nocall_haploid(self):
        nucleotide_genotypes = [('D', 'I'), ('I', 'D')]
        result = convert_indel_genotype_to_vcf(nucleotide_genotypes, DummyRecord(), True, 1)
        self.assertEqual(result, ".")

    def test_consistent_DI_DI_returns_0_1_deletion_true(self):
        nucleotide_genotypes = [('D', 'I'), ('D', 'I')]
        result = convert_indel_genotype_to_vcf(nucleotide_genotypes, DummyRecord(), True, 2)
        self.assertEqual(result, "0/1")

    def test_consistent_DI_DI_returns_0_1_deletion_false(self):
        nucleotide_genotypes = [('D', 'I'), ('D', 'I')]
        result = convert_indel_genotype_to_vcf(nucleotide_genotypes, DummyRecord(), False, 2)
        self.assertEqual(result, "0/1")

if __name__ == "__main__":
    unittest.main()