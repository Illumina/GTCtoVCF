import os
import shutil
import sys
import tempfile
import unittest
from logging import Logger
from subprocess import call

from gtc_to_vcf.BPMRecord import BPMRecord
from gtc_to_vcf.BPMRecord import IndelSourceSequence, determine_left_shift
from gtc_to_vcf.GenotypeFormat import RecordCombiner
from gtc_to_vcf.IlluminaBeadArrayFiles import RefStrand

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, os.pardir))


class Regression(unittest.TestCase):
    @staticmethod
    def get_script():
        return os.path.join(SCRIPT_DIR, "..", "gtc_to_vcf.py")

    @staticmethod
    def get_genome():
        return os.path.join(SCRIPT_DIR, "data", "genome.fa")

    def run_regression(self, command, output, expected_output):
        ret_code = call(command)
        self.assertEqual(ret_code, 0)
        self.compare_vcf(output, expected_output)

    def compare_vcf(self, output, expected_output):
        self.assertTrue(os.path.isfile(output))
        with open(output) as file1, open(expected_output) as file2:
            for (line1, line2) in zip(file1, file2):
                if line1.startswith("##source") and line2.startswith("##source"):
                    continue
                if line1.startswith("##reference") and line2.startswith("##reference"):
                    continue
                self.assertEqual(line1, line2)


class RegressionBPM(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.bpm"), "--output-vcf-path", output_vcf,
                   "--skip-indels", "--disable-genome-cache"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionBPM", "output", "output.vcf"))
        os.remove(output_vcf)


class RegressionCSV(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--output-vcf-path", output_vcf,
                   "--disable-genome-cache"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionCSV", "output", "output.vcf"))
        os.remove(output_vcf)


class RegressionLociFilter(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--filter-loci",
                   os.path.join(SCRIPT_DIR, "data", "RegressionLociFilter", "input", "loci_file.txt"),
                   "--output-vcf-path", output_vcf, "--disable-genome-cache"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionLociFilter", "output", "output.vcf"))
        os.remove(output_vcf)


class RegressionGTC(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--gtc-paths",
                   os.path.join(SCRIPT_DIR, "data", "RegressionGTC", "input", "201142750001_R01C01.gtc"),
                   "--output-vcf-path", output_vcf, "--disable-genome-cache"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionGTC", "output", "output.vcf"))
        os.remove(output_vcf)


class RegressionGTCDirectory(Regression):
    def test(self):
        output_vcf_dir = tempfile.mkdtemp()
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--gtc-paths",
                   os.path.join(SCRIPT_DIR, "data", "RegressionGTCDirectory", "input"), "--output-vcf-path",
                   output_vcf_dir, "--disable-genome-cache"]
        self.run_regression(command, os.path.join(output_vcf_dir, "201142750001_R01C01.vcf"),
                            os.path.join(SCRIPT_DIR, "data", "RegressionGTCDirectory", "output",
                                         "201142750001_R01C01.vcf"))
        self.run_regression(command, os.path.join(output_vcf_dir, "201142750001_R01C02.vcf"),
                            os.path.join(SCRIPT_DIR, "data", "RegressionGTCDirectory", "output",
                                         "201142750001_R01C02.vcf"))
        shutil.rmtree(output_vcf_dir)


class RegressionAux(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--output-vcf-path", output_vcf,
                   "--auxiliary-loci", os.path.join(SCRIPT_DIR, "data", "RegressionAux", "input", "auxiliary.vcf"),
                   "--disable-genome-cache"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionAux", "output", "output.vcf"))
        os.remove(output_vcf)


class RegressionIncludeAttributes(Regression):
    def test(self):
        output_vcf = tempfile.mktemp(suffix=".vcf")
        command = [sys.executable, self.get_script(), "--genome-fasta-file", self.get_genome(), "--manifest-file",
                   os.path.join(SCRIPT_DIR, "data", "small_manifest.csv"), "--gtc-paths",
                   os.path.join(SCRIPT_DIR, "data", "RegressionIncludeAttributes", "input", "201142750001_R01C01.gtc"),
                   "--output-vcf-path", output_vcf, "--disable-genome-cache", "--include-attributes", "GT", "GQ", "BAF",
                   "LRR"]
        self.run_regression(command, output_vcf,
                            os.path.join(SCRIPT_DIR, "data", "RegressionIncludeAttributes", "output", "output.vcf"))
        os.remove(output_vcf)


class TestSourceSequence(unittest.TestCase):
    def test_source_sequence_split(self):
        (five_prime, indel, three_prime) = IndelSourceSequence.split_source_sequence("ACGT[-/ATAT]GGTA")
        self.assertEqual(five_prime, "ACGT")
        self.assertEqual(indel, "ATAT")
        self.assertEqual(three_prime, "GGTA")

        (five_prime, indel, three_prime) = IndelSourceSequence.split_source_sequence("[-/ATAT]")
        self.assertEqual(five_prime, "")
        self.assertEqual(indel, "ATAT")
        self.assertEqual(three_prime, "")

        self.assertEqual(determine_left_shift("AACAA", "A", "GGG"), ("AAC", "AAGGG"))
        self.assertEqual(determine_left_shift("AACAA", "AAA", "GGG"), ("AAC", "AAGGG"))
        self.assertEqual(determine_left_shift("AACAA", "AAA", "AGGG"), ("AAC", "AAAGGG"))


class TestCombinedGenotypes(unittest.TestCase):
    def check_genotype(self, data, expected_genotype):
        records = [datum[0] for datum in data]
        for idx in range(len(records)):
            records[idx].index_num = idx
        genotypes = [datum[1] for datum in data]
        combiner = RecordCombiner(records, genotypes, "")
        self.assertEqual(sorted(combiner.combine_genotypes()), sorted(expected_genotype))

    # noinspection PyTypeChecker
    def test_genotype_combination(self):
        logger = Logger("test_genotype_combinations")
        # Inf II [A/C] -> AC
        # AC
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 0, logger), 2)]
        self.check_genotype(data, ('A', 'C'))

        # Inf II [A/C] -> AC
        # Inf II [A/C] -> CC
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 0, logger), 2),
                (BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 0, logger), 3)]
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/C] -> AA
        # AA
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 1)]
        self.check_genotype(data, ('A', 'A'))

        # Inf II [A/C] -> AC
        # Inf II [A/C] -> NC
        # AC
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 0)]
        self.check_genotype(data, ('A', 'C'))

        # Inf II [A/C] -> AC
        # Inf I [A/G] -> NC
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 0)]
        self.check_genotype(data, ('-', '-'))

        # Inf I [A/C] -> AA
        # Inf I [A/G] -> AG
        # AG
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 1),
                (BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 2)]
        self.check_genotype(data, ('A', 'G'))

        # Inf I [A/C] -> AC
        # Inf I [A/G] -> AG
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/C]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 2)]
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> AA
        # Inf I [T/G] -> NC
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 1),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 0)]
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> GG
        # AG
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 3)]
        self.check_genotype(data, ('A', 'G'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> GG
        # AG
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 3)]
        self.check_genotype(data, ('A', 'G'))

        # Inf II [A/G] -> AA
        # Inf I [T/G] -> TT
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 1),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 1)]
        self.check_genotype(data, ('-', '-'))

        # Inf II [A/G] -> GG
        # Inf I [T/G] -> GG
        # GG
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 3),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 3)]
        self.check_genotype(data, ('G', 'G'))

        # Inf II [A/G] -> AG
        # Inf I [T/G] -> TG
        # TG
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 2),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 2)]
        self.check_genotype(data, ('T', 'G'))

        # Inf II [A/G] -> AA
        # Inf I [T/G] -> TT
        # --
        data = [(BPMRecord("", None, None, "", 0, "[A/G]", RefStrand.Plus, 0, None, None, None, None, 1, logger), 1),
                (BPMRecord("", None, None, "", 0, "[T/G]", RefStrand.Plus, 1, None, None, None, None, 1, logger), 1)]
        self.check_genotype(data, ('-', '-'))


def main():
    reference_file = os.path.join(SCRIPT_DIR, "data", "genome.fa")
    if not os.path.isfile(reference_file):
        download_command = ["bash", os.path.join(SCRIPT_DIR, os.path.pardir, "scripts", "download_reference.sh"),
                            reference_file]
        call(download_command)
    unittest.main()


if __name__ == "__main__":
    main()
