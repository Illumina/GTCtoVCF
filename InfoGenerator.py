class LocusInfo(object):
    """
    Stores locus info for a single entry.
    Currently just a place-holder
    """

    def __init__(self):
        """
        Does nothing
        """
        pass

    def update(self, vcf_record):
        """
        Update VCF record with finalized information

        Currently does nothing
        """
        pass

    def add_sample(self, bpm_records, gtc):
        """
        Update locus info with information from sample

        Currently does nothing
        """
        pass


class LocusInfoFactory(object):
    """
    Class to create locus info objects. Also will return the neccessary
    VCF header information related to the info fields. Configuring this
    class will allow specification of different Locus Info fields
    """

    def __init__(self):
        pass

    def get_vcf_headers(self):
        """
        Returns list of VCF info header lines
        """
        return []

    def create_locus_info(self):
        """
        Return new LocusInfo object
        """
        return LocusInfo()
