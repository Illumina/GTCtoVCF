from IlluminaBeadArrayFiles import BeadPoolManifest, RefStrand
from BPMRecord import BPMRecord, reverse_complement

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
            yield BPMRecord(bpm.names[idx], bpm.addresses[idx], bpm.chroms[idx], bpm.map_infos[idx], bpm.snps[idx], bpm.ref_strands[idx], bpm.assay_types[idx], None, idx)


def get_indel_sequence(source_sequence):
    """
    Extract indel sequence from a source sequence string

    Args:
        source_sequence (string): Source sequence string

    Returns:
        string: The sequence of the indel
    """
    left_position = source_sequence.find("/")
    right_position = source_sequence.find("]")
    assert source_sequence[left_position - 1] == "-"
    return source_sequence[(left_position + 1):right_position]


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
                                  "chr", "mapinfo", "refstrand", "sourceseq", "snp", "addressb_id")
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
                (source_strand, ilmn_strand, name, chrom, map_info, ref_strand, source_seq, snp, addressb_id) = [
                    bits[required_column2idx[column]] for column in self._required_columns]

                #
                # determine the source snp on the design strand
                #
                source_strand = source_strand[0].upper()
                ilmn_strand = ilmn_strand[0].upper()

                if "D" in snp:
                    if source_strand == "U" or ilmn_strand == "U":
                        logger.warn("Skipping indel " + name +
                                    " with unknown ILMN or source strand")
                        continue

                    indel_sequence = get_indel_sequence(source_seq)
                    if source_strand == "P" or source_strand == "M":
                        assert ilmn_strand == "P" or ilmn_strand == "M"
                    else:
                        assert ilmn_strand == "T" or ilmn_strand == "B"

                    if source_strand != ilmn_strand:
                        indel_sequence = reverse_complement(indel_sequence)
                else:
                    indel_sequence = None

                assay_type = 0 if addressb_id == "" else 1
                yield BPMRecord(name, 0, chrom, map_info, snp, RefStrand.from_string(ref_strand), assay_type, indel_sequence, idx)
