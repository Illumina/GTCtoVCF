import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from GenotypeGenerator import RecordCombiner
from BPMRecord import BPMRecord
from IlluminaBeadArrayFiles import RefStrand
from BPMRecord import split_source_sequence

class TestSourceSequence(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_source_sequence_split(self):
        (five_prime, indel, three_prime) = split_source_sequence("ACGT[-/ATAT]GGTA")
        self.assertEqual(five_prime, "ACGT")
        self.assertEqual(indel, "ATAT")
        self.assertEqual(three_prime, "GGTA")

        (five_prime, indel, three_prime) = split_source_sequence("[-/ATAT]")
        self.assertEqual(five_prime, "")
        self.assertEqual(indel, "ATAT")
        self.assertEqual(three_prime, "")

class TestCombinedGenotypes(unittest.TestCase):
    def setUp(self):
		# read in multiple SNP bpm
		# create record objects
		# look for the combined name in the record
		# check that the genotypes agree
        pass

    def check_genotype(self, data, expected_genotype):
        records = [datum[0] for datum in data]
        for idx in xrange(len(records)):
            records[idx].index_num = idx
        genotypes = [datum[1] for datum in data]
        combiner = RecordCombiner(records, genotypes, "")
        self.assertEqual(combiner.combine_genotypes(), expected_genotype)

    def test_genotype_combination(self):
        # Inf II [A/C] -> AC
        # AC
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 0), 2))
        self.check_genotype(data, ('A', 'C'))

        # Inf II [A/C] -> AC
        # Inf II [A/C] -> CC
        # --
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 0), 2))
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 0), 3))
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/C] -> AA
        # AA
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 1), 1))
        self.check_genotype(data, ('A', 'A'))

        # Inf II [A/C] -> AC
        # Inf II [A/C] -> NC
        # AC
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 1), 0))
        self.check_genotype(data, ('A', 'C'))

        # Inf II [A/C] -> AC
        # Inf I [A/G] -> NC
        # --
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 0, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 1, "", 1), 0))
        self.check_genotype(data, ('-', '-'))

        # Inf I [A/C] -> AA
        # Inf I [A/G] -> AG
        # AG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 1, "", 1), 1))
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 1, "", 1), 2))
        self.check_genotype(data, ('A', 'G'))

        # Inf I [A/C] -> AC
        # Inf I [A/G] -> AG
        # --
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/C]", RefStrand.Plus, 1, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 1, "", 1), 2))
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> AA
        # Inf I [T/G] -> NC
        # --
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 1))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 0))
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> GG
        # AG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 3))
        self.check_genotype(data, ('A', 'G'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> GG
        # AG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 3))
        self.check_genotype(data, ('A', 'G'))

        # Inf II [A/G] -> AA
        # Inf I [T/G] -> TT
        # --
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 1))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 1))
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> GG
        # Inf I [T/G] -> GG
        # GG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 3))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 3))
        self.check_genotype(data, ('G', 'G'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> TG
        # TG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 2))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 2))
        self.check_genotype(data, ('T', 'G'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> TG
        # TG
        data = []
        data.append((BPMRecord("", None, "", 0, "[A/G]", RefStrand.Plus, 0, "", 1), 1))
        data.append((BPMRecord("", None, "", 0, "[T/G]", RefStrand.Plus, 1, "", 1), 1))
        self.check_genotype(data, ('-', '-'))

if __name__ == "__main__":
    unittest.main()

