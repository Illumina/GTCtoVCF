from IlluminaBeadArrayFiles import BeadPoolManifest, RefStrand
from BPMRecord import BPMRecord, reverse_complement, IndelSourceSequence

class BPMReader(object):
    """Get records from a BPM manifest file"""

    def __init__(self, bpm_file):
        """
        Initialize a BPM reader with a file path

        Args:
            bpm_file (string): Path to the BPM manifest

        Returns:
            BeadPoolReader
        """
        self.bpm = BeadPoolManifest(bpm_file)

    def get_bpm_records(self):
        """
        Get BPM records from the reader

        Args:
            None

        Yields:
            BPMRecord: Next BPMRecord in the file
        """
        bpm = self.bpm
        for idx in xrange(len(bpm.addresses)):
            yield BPMRecord(bpm.names[idx], bpm.addresses[idx], None, bpm.chroms[idx], bpm.map_infos[idx], bpm.snps[idx], bpm.ref_strands[idx], bpm.assay_types[idx], None, idx)


class CSVManifestReader(object):
    """Get records from a CSV manifest"""

    def __init__(self, csv_file, logger):
        """
        Initialize a manifest reader from a CSV file

        Args:
            csv_file (string): Path to the CSV manifest
            logger (logging.Logger): Logger for reporting logging information

        Returns:
            CSVManifestReader
        """
        self._csv_file = csv_file
        self._required_columns = ("sourcestrand", "ilmnstrand", "name",
                                  "chr", "mapinfo", "refstrand", "sourceseq", "snp", "addressb_id", "allelea_probeseq")
        self._logger = logger

    def get_bpm_records(self):
        """
        Get BPM records from the reader

        Args:
            None

        Yields:
            BPMRecord: Next BPMRecord in the file

        Raises:
            Exception - Manifest is missing required column
        """
        in_data = False
        idx = -1
        for line in open(self._csv_file):
            if line.startswith("IlmnID,"):
                in_data = True
                header = line.rstrip().lower().split(",")
                required_column2idx = {}
                for required_column in self._required_columns:
                    try:
                        required_column2idx[required_column] = header.index(
                            required_column)
                    except:
                        raise Exception(
                            "Manifest is missing required column " + required_column)
                continue
            if line.startswith("[Controls]"):
                in_data = False
                continue

            if in_data:
                idx += 1
                bits = line.rstrip().split(",")
                (source_strand, ilmn_strand, name, chrom, map_info, ref_strand, source_seq, snp, addressb_id, probe_a) = [
                    bits[required_column2idx[column]] for column in self._required_columns]

                #
                # determine the source snp on the design strand
                #
                source_strand = source_strand[0].upper()
                ilmn_strand = ilmn_strand[0].upper()

                if "D" in snp:
                    indel_source_sequence = IndelSourceSequence(source_seq, source_strand, ilmn_strand)
                else:
                    indel_source_sequence = None

                assay_type = 0 if addressb_id == "" else 1
                yield BPMRecord(name, 0, probe_a, chrom, map_info, snp, RefStrand.from_string(ref_strand), assay_type, indel_source_sequence, idx)
