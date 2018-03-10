from IlluminaBeadArrayFiles import BeadPoolManifest, RefStrand
from BPMRecord import BPMRecord, IndelSourceSequence

class BPMReader(object):
    """
    Get records from a BPM manifest file

    Attributes:
        source_file (string) : Source file used to create reader
    """

    def __init__(self, bpm_file):
        """
        Initialize a BPM reader with a file path

        Args:
            bpm_file (string): Path to the BPM manifest

        Returns:
            BeadPoolReader
        """
        self.source_file = bpm_file
        self._bpm = BeadPoolManifest(bpm_file)

    def get_bpm_records(self):
        """
        Get BPM records from the reader

        Args:
            None

        Yields:
            BPMRecord: Next BPMRecord in the file
        """
        bpm = self._bpm
        for idx in xrange(len(bpm.addresses)):
            yield BPMRecord(bpm.names[idx], bpm.addresses[idx], None, bpm.chroms[idx], bpm.map_infos[idx], bpm.snps[idx], bpm.ref_strands[idx], bpm.assay_types[idx], None, None, None, None, idx)

class ManifestFilter(object):
    """
    Filter entries from another manifest reader

    Attributes:
        source_file (string) : The file used a source for records
    """

    def __init__(self, manifest_reader, loci_to_filter, skip_indels, logger):
        """
        Return a new ManifestFilter. Will skip records as specified in constructor
        as well as records with chromosome or mapping of zero.

        Args:
            manifest_reader (CSVManifestReader/BPMReader) : The source of BPM records
            loci_to_filter (set(string)) : A set of record names to skip, may be None
            skip_indels (bool) : Skip indels
            logger (Logger) : logger

        Returns:
            ManifestFilter
        """
        self._manifest_reader = manifest_reader
        self._loci_to_filter = loci_to_filter
        self._skip_indels = skip_indels
        self._logger = logger
        self.source_file = manifest_reader.source_file

    def get_bpm_records(self):
        for record in self._manifest_reader.get_bpm_records():
            if record.chromosome == "0" or record.pos == 0:
                continue
            if self._loci_to_filter and record.name in self._loci_to_filter:
                continue
            if record.is_indel() and self._skip_indels:
                self._logger.warning("Skipping indel " + record.name)
                continue
            yield record

class CSVManifestReader(object):
    """
    Get records from a CSV manifest

    Attributes:
        source_file (string) : Source file used to create reader
    """

    def __init__(self, csv_file, genome_reader, logger):
        """
        Initialize a manifest reader from a CSV file

        Args:
            csv_file (string): Path to the CSV manifest
            genome_reader (ReferenceGenome,CachedReferenceGenome)
            logger (logging.Logger): Logger for reporting logging information

        Returns:
            CSVManifestReader
        """
        self.source_file = csv_file
        self._required_columns = ("sourcestrand", "ilmnstrand", "name",
                                  "chr", "mapinfo", "refstrand", "sourceseq", "snp", "addressb_id", "allelea_probeseq")
        self._genome_reader = genome_reader
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
        for line in open(self.source_file):
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

                if "D" in snp:
                    indel_source_sequence = IndelSourceSequence(source_seq)
                else:
                    indel_source_sequence = None

                assay_type = 0 if addressb_id == "" else 1
                try:
                    yield BPMRecord(name, 0, probe_a, chrom, map_info, snp, RefStrand.from_string(ref_strand), assay_type, indel_source_sequence, source_strand, ilmn_strand, self._genome_reader, idx)
                except Exception as error:
                    self._logger.warn("Failed to process entry for record %s: %s", name, str(error))
