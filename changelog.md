# 1.2.1 (2020/03/28)
* Fix for large manifests with norm id byte-wrapping issue (e.g. Omni5)
* Fixes for aux VCF MNV handling
* Updates GenCallScore to phred-scaled int
* Updates to download_reference.sh for new genome build location and to work on OSX

# 1.2.0 (2019/05/24)
* Fixes for python3 compatibility
* Use chromosome ordering from reference
* Use semi-colon as delimiter within identifier field
* Fix for download_reference.sh to support more UNIX-y platforms
* Various small typos
* Support for outputting BAF and LRR

# 1.1.1 (2018/03/27)
* Output warning if any mismatch to reference genome in indel source sequence
* Count length of indel when scoring deletion reference allele

# 1.1.0 (2018/03/20)
* Enable caching of genomic reference data
* Enable processing of multiple GTC files or directories for improved performance
* Address issue with incorrect assignment of reference allele for indels
* Added script to build reference FASTA files for human

# 1.0.1 (2017/10/27)
* Eliminated extraneous header in column

# 1.0.0 (2017/10/15)
* Initial release
