from vcf.parser import _Format
import numpy as np


class LogRRatioFormat(object):
    """
    Generate log R ratio format information for VCF
    """
    def __init__(self, logger, log_r_ratio):
        self._log_r_ratio = log_r_ratio
        self._logger = logger

    @staticmethod
    def get_id():
        return "LRR"

    @staticmethod
    def get_description():
        return "Log R Ratio"

    @staticmethod
    def get_format_obj():
        # arguments should be: ['id', 'num', 'type', 'desc']
        return _Format(LogRRatioFormat.get_id(), 1, "Float", LogRRatioFormat.get_description())

    def generate_sample_format_info(self, bpm_records, vcf_record, sample_name):
        """
        Get the sample log R ratio

        Args:

        Returns:
            float: log R ratio
        """
        # if we have more than 1 BPMRecord, need to merge
        log_r_ratio_list = []
        if len(bpm_records) > 1:
            # if we have more than 1 alt allele, return missing
            if len(vcf_record.ALT) > 1:
                return "."
            for i in range(len(bpm_records)):
                idx = bpm_records[i].index_num
                log_r_ratio = self._log_r_ratio[idx]
                log_r_ratio_list.append(log_r_ratio)

                # snp = bpm_records[i].snp
                # strand = bpm_records[i].ref_strand
                # ref_allele = vcf_record.REF

                # chrom = bpm_records[i].chromosome
                # start_pos = bpm_records[i].pos
                # alt_allele = vcf_record.ALT
                # print("chr: {} pos: {} ref: {} alt: {} snp: {} strand: {} LRR: {}".format(chrom, start_pos, ref_allele, alt_allele, snp, strand, log_r_ratio))

        # otherwise single BPMRecord
        else:
            # if we have a multi-allelic site, return missing
            if len(vcf_record.ALT) > 1:
                return "."
            idx = bpm_records[0].index_num
            log_r_ratio_list.append(self._log_r_ratio[idx])

        # nanmedian ignores NaN values
        final_log_r_ratio = np.nanmedian(log_r_ratio_list)
        if np.isnan(final_log_r_ratio):
            return "."
        else:
            return final_log_r_ratio
