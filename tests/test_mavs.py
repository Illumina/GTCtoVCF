import os
import sys
import unittest
from logging import Logger

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, os.pardir))

from GenotypeFormat import RecordCombiner
from BPMRecord import BPMRecord
from IlluminaBeadArrayFiles import RefStrand

class TestCombinedGenotypesMavs(unittest.TestCase):
    def check_genotype(self, data, expected_genotype):
        records = [datum[0] for datum in data]
        for idx in range(len(records)):
            records[idx].index_num = idx
        genotypes = [datum[1] for datum in data]
        combiner = RecordCombiner(records, genotypes, "")
        self.assertEqual(sorted(combiner.combine_genotypes()), sorted(expected_genotype))

    def case1(self):
        # MAV = T -> C,G
        # Inf II [T/C] -> CC
        # Inf I [T/C] -> CC
        # Inf I [T/G] -> NC
        # Combined Genotype -> NC
        logger = Logger("test_genotype_combinations")
        data = []
        # Note the last field is the genotype where (0 = NC, 1 = AA, 2 = AB, 3 = BB)
        # Note the field after the RefStrand is the assay_type where (0 = InfII, 1 = InfI)
        data.append((BPMRecord("ilmnseq_2:234627937", None, None, "", 0, "[A/G]", RefStrand.Minus, 0, None, None, None, None, 0, logger), 3))
        data.append((BPMRecord("ilmnseq_rs2011404.1_F2BT", None, None, "", 0, "[A/G]", RefStrand.Minus, 1, None, None, None, None, 0, logger), 3))
        data.append((BPMRecord("ilmnseq_rs2011404.2_F2BT", None, None, "", 0, "[A/C]", RefStrand.Minus, 1, None, None, None, None, 0, logger), 0))

        self.check_genotype(data, ('-', '-'))

if __name__ == "__main__":
    unittest.main()