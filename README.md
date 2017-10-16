# GTC to VCF converter

## Usage
```bash
usage: gtc_to_vcf.py [-h] [--gtc-file GTC_FILE] --manifest-file MANIFEST_FILE
                     --genome-fasta-file GENOME_FASTA_FILE
                     [--output-vcf-file OUTPUT_VCF_FILE] [--skip-indels]
                     [--log-file LOG_FILE] [--expand-identifiers]
                     [--unsquash-duplicates] [--auxiliary-loci AUXILIARY_LOCI]
                     [--version]

Convert GTC file to VCF format

optional arguments:
  -h, --help            show this help message and exit
  --gtc-file GTC_FILE   GTC file
  --manifest-file MANIFEST_FILE
                        Bead pool manifest for product (*.csv or *.bpm)
  --genome-fasta-file GENOME_FASTA_FILE
                        Reference genome in fasta format
  --output-vcf-file OUTPUT_VCF_FILE
                        Name of VCF file to be created (default is output.vcf)
  --skip-indels         Skip processing of indels (default is False)
  --log-file LOG_FILE   File to write logging information
  --expand-identifiers  For VCF entries with multiple corresponding manifest
                        entries, list all manifest identifiers in VCF ID field
  --unsquash-duplicates
                        Generate unique VCF records for duplicate assays
  --auxiliary-loci AUXILIARY_LOCI
                        VCF file with auxiliary definitions of loci
  --version             show program's version number and exit
```
## Input details
### Manifests
The supplied manifest file may either be in CSV or BPM format; however, a CSV format manifest is required to generate indel records in the output VCF. When running the converter with a BPM manifest, indel processing must be explicitly disabled with the "--skip-indels" option. In either case, the manifest must provide RefStrand annotations. The GTC to VCF converter depends on the presence of accurate mapping information within the manifest, and may produce inaccurate results if the mapping information is incorrect. Mapping information should follow the implicit dbSNP standard, where
* Positions are reported with 1-based indexing
* Positions in the PAR are reported with mapping position to the X chromosome
* For an insertion relative to the reference, the position of the base immediately 5' to the insertion (on the plus strand) is given
* For a deletion relative to the reference, the position of most 5' deleted based (on the plus strand) is given

Any standard product manifest provided by Illumina will already follow these conventions. 
### Reference genome
The contig identifiers in the provided genome FASTA file must match exactly the chromosome identifiers specified in the provided manifest. For a standard human product manifest, this means that the contig headers should read ">1" rather than ">chr1". For compatibility with BaseSpace Variant Interpreter (https://www.illumina.com/informatics/research/biological-data-interpretation/variant-interpreter.html), the specified path of the reference must contain either contain the string GrCh37 or GrCh38, for build 37 and 38, respectively. Suitable whole genome FASTA files are provided in the following Illumina iGenomes builds: ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz (for build 37) and ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz (for build 38).

### Auxiliary loci
Certain classes of variant types (such as multi-nucleotide variants) are not currently supported in the upstream analysis software that produces GTC files. However, it is possible to query this type of variant by creating a SNP design that differentiates the specific multi-nucleotide alleles of interest. For example, if the true source sequence is

ATGC[AT/CG]GTAA

This assay could be designed as a SNP assay with the following source sequence

ATGC[A/C]NNNN

The GTC converter tool provides an option to supply a list of auxiliary records (in VCF format) to restore the true alleles for these cases in the output VCF. There are several restrictions around this function
* The auxiliary definition must be bi-allelic
* The auxiliary definition must be a multi-nucleotide variant
* There must not be multiple array assays (e.g., duplicates) for the locus 

## Output description
The VCF file output follows VCF4.1 format (https://samtools.github.io/hts-specs/VCFv4.1.pdf). Some additional details on output formatting:
* Genotypes are adjusted to reflect the sample ploidy. Calls are haploid for loci on Y, MT and non-PAR chromosome X for males.
* Multiple SNPs in the input manifest which are mapped to the same chromosomal coordinate (e.g. tri-allelic loci or duplicated sites) are collapsed into one VCF entry and a combined genotype generated. To produce the combined genotype, the set of all possible genotypes is enumerated based on the queried alleles. Genotypes which are not possible based on called alleles and assay design limitations (e.g. infiniumII designs cannot distinguish between A/T and C/G calls) are filtered. If only one consistent genotype remains after the filtering process, then the site is assigned this genotype  Otherwise, the genotype is ambiguous (more than 1) or inconsistent (less than 1) and a no-call is returned.  Please see 'test_class.py' in the 'tests' folder for unit tests and examples demonstrating  the genotype merging process.

## Requirements
The following python packages are required. Note that the pysam package is currently not supported on the Windows platform.

Requirement | Version
--- | ---
pysam | 0.9.0
numpy | 1.11.2
pyvcf | 0.6.8

An easy way to obtain a version of python with these dependencies available is to use "Miniconda". Miniconda is a minimal python installation along with a package manager ("conda") that can be used to install additional packages. First, obtain the 64-bit python 2.7 installer from https://conda.io/miniconda.html and install
```bash
bash Miniconda2-latest-Linux-x86_64.sh 
```
This will run the installer and ask you where you would like to install Miniconda. Then, install numpy, pyvcf, and pysam
```bash
conda install -c miniconda numpy=1.11.2
conda install -c bioconda pyvcf=0.6.8
conda install -c bioconda pysam=0.9.0
```
where conda is a the package manager binary located in the installation location specified in the first step. 

## License

Copyright 2017 Illumina

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
