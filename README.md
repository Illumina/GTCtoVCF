# GTC to VCF converter

## Semi-archived state
As of May 2023, it is now recommended to use the fully-supported Illumina solution:
[DRAGEN Array](https://support.illumina.com/array/array_software/dragen-array-secondary-analysis.html) for GTC to VCF conversion. That tool is the most up-to-date with bug fixes and improvements.
This repository will remain open, but may be fully archived in the near future.

## Requirements

**GTCtoVCF is currently only known to work with python2.7.**
The following packages are also required.
Note that the pysam package is currently not supported on the Windows platform.

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


## Usage

```none
usage: gtc_to_vcf.py [-h] [--gtc-paths GTC_PATHS [GTC_PATHS ...]]
                     --manifest-file MANIFEST_FILE --genome-fasta-file
                     GENOME_FASTA_FILE [--output-vcf-path OUTPUT_VCF_PATH]
                     [--skip-indels] [--log-file LOG_FILE]
                     [--expand-identifiers] [--unsquash-duplicates]
                     [--auxiliary-loci AUXILIARY_LOCI]
                     [--filter-loci FILTER_LOCI] [--disable-genome-cache]
                     [--include-attributes [{GT,GQ,BAF,LRR} [{GT,GQ,BAF,LRR} ...]]]
                     [--version]

Convert GTC file to VCF format

optional arguments:
  -h, --help            show this help message and exit
  --gtc-paths GTC_PATHS [GTC_PATHS ...]
                        One or more GTC files or directories to process (optional)
  --manifest-file MANIFEST_FILE
                        Bead pool manifest for product (*.csv or *.bpm)
  --genome-fasta-file GENOME_FASTA_FILE
                        Reference genome in fasta format
  --output-vcf-path OUTPUT_VCF_PATH
                        Path for generation of VCF output (default is
                        output.vcf)
  --skip-indels         Skip processing of indels (default is False)
  --log-file LOG_FILE   File to write logging information (optional)
  --expand-identifiers  For VCF entries with multiple corresponding manifest
                        entries, list all manifest identifiers in VCF ID field
  --unsquash-duplicates
                        Generate unique VCF records for duplicate assays
  --auxiliary-loci AUXILIARY_LOCI
                        VCF file with auxiliary definitions of loci (optional)
  --filter-loci FILTER_LOCI
                        File containing list of loci names to filter from
                        input manifest (optional)
  --disable-genome-cache
                        Disable caching of genome reference data
  --include-attributes [{GT,GQ,BAF,LRR} [{GT,GQ,BAF,LRR} ...]]
                        Additional attributes to include in VCF FORMAT output
                        (optional) (default is GT GQ)
  --version             show program's version number and exit


```
## Input details
### Input and output files
Input GTC files are specified with the --gtc-paths option. One or more paths may be specified with this option. If a given path corresponds to a directory, the script will identify all GTC files within that directory for processing. A combination of file and directory paths may be specified. There is a performance benefit to analyzing multiple GTC files in a single invocation of the program, as it minimizes the IO overhead of reading manifest and genome reference data to memory. When multiple input GTC files are specified, a single VCF file is produced for each input GTC file, as opposed to a single multi-sample VCF file. When no GTC file is provided, the program will still produce an output VCF file without sample genotyping information.

The --output-vcf-path option may either be a file or directory. If the argument is a directory, the program will automatically determine an appropriate name for the VCF output file created within that directory. The behavior is summarized in the following table

| # GTC files | --vcf-output-path | Behavior |
| - | - | - |
| 0 | directory  | VCF filename determined from manifest filename  |
| 1 | directory  | VCF filename determined from input GTC file  |
| 2+ | directory  | VCF filenames determined from input GTC files  |
| 0 | file  | VCF filename determined from --vcf-output-path argument  |
| 1 | file  | VCF filename determined from --vcf-output-path argument  |
| 2+ | file  | Error |

### Manifests
The supplied manifest file may either be in CSV or BPM format; however, a CSV format manifest is required to generate indel records in the output VCF. When running the converter with a BPM manifest, indel processing must be explicitly disabled with the "--skip-indels" option. In either case, the manifest must provide RefStrand annotations. The GTC to VCF converter depends on the presence of accurate mapping information within the manifest, and may produce inaccurate results if the mapping information is incorrect. Mapping information should follow the implicit dbSNP standard, where
* Positions are reported with 1-based indexing.
* Positions in the PAR are reported with mapping position to the X chromosome.
* For an insertion relative to the reference, the position of the base immediately 5' to the insertion (on the plus strand) is given.
* For a deletion relative to the reference, the position of the most 5' deleted based (on the plus strand) is given.

Any standard product manifest provided by Illumina will already follow these conventions.
### Reference genome
The contig identifiers in the provided genome FASTA file must match exactly the chromosome identifiers specified in the provided manifest. For a standard human product manifest, this means that the contig headers should read ">1" rather than ">chr1". For compatibility with BaseSpace Variant Interpreter (https://www.illumina.com/informatics/research/biological-data-interpretation/variant-interpreter.html), the specified path of the reference must contain either contain the string GrCh37 or GrCh38, for build 37 and 38, respectively. Suitable whole genome FASTA files can be built with the download_reference.sh script located within the scripts directory. This bash script is dependent on samtools (http://www.htslib.org/download/). If running on OSX, you may need to install coreutils (e.g. `brew install coreutils`).

### Squashing duplicates
In the manifest, there can be cases where the same variant is probed by multiple different assays. These assays may be the same design or alternate designs for the same locus. In the default mode of operation, these duplicates will be "squashed" into a single record in the VCF. The method used to incorporate information across multiple assays is under the latter "Output description" heading. When the "--unsquash-duplicates" option is provided, this "squashing" behavior is disabled, and each duplicate assay will be reported in a separate entry in the VCF file. This option is helpful when you are interested in investigating or validating the performance of individual assays, rather than trying to generate genotypes for specific variants. Note that if a locus has more than two alleles and is also queried with duplicated designs, the duplicates will not be unsquashed.

### Genome cache
By default, the entire reference genome will be read into memory. Generally, this will be more efficient than reading data from the indexed reference on disk at the expense of greater memory utilization. For situations in which the genome caching is not desirable (low memory availability or a small input manifest), it is possible to disable this default behavior with the "--disable-genome-cache" option.

### Auxiliary loci
Certain classes of variant types (such as multi-nucleotide variants) are not currently supported in the upstream analysis software that produces GTC files. However, it is possible to query this type of variant by creating a SNP design that differentiates the specific multi-nucleotide alleles of interest. For example, if the true source sequence is

ATGC[AT/CG]GTAA

This assay could be designed as a SNP assay with the following source sequence

ATGC[A/C]NNNN

The GTC converter tool provides an option to supply a list of auxiliary records (in VCF format) to restore the true alleles for these cases in the output VCF. There are several restrictions around this function
* The auxiliary definition must be bi-allelic.
* The auxiliary definition must be a multi-nucleotide variant.
* There must not be multiple array assays (e.g., duplicates) for the locus.

### Include Attributes
Default attributes to be included in the VCF are "GT" (genotype) and "GQ" (genotype score).
With this option, you can compose which fields from the GTC are output to the VCF.
The additional fields available (from GTC version 5) are "BAF" (B Allele Frequency) and "LRR" (Log R Ratio).

## Output description
The VCF file output follows VCF4.1 format (https://samtools.github.io/hts-specs/VCFv4.1.pdf). Some additional details on output formatting:
* Genotypes are adjusted to reflect the sample ploidy. Calls are haploid for loci on Y, MT and non-PAR chromosome X for males.
* Multiple SNPs in the input manifest which are mapped to the same chromosomal coordinate (e.g. tri-allelic loci or duplicated sites) are collapsed into one VCF entry and a combined genotype generated. To produce the combined genotype, the set of all possible genotypes is enumerated based on the queried alleles. Genotypes which are not possible based on called alleles and assay design limitations (e.g. infiniumII designs cannot distinguish between A/T and C/G calls) are filtered. If only one consistent genotype remains after the filtering process, then the site is assigned this genotype  Otherwise, the genotype is ambiguous (more than 1) or inconsistent (less than 1) and a no-call is returned.  Please see 'test_class.py' in the 'tests' folder for unit tests and examples demonstrating  the genotype merging process.

## Docker

### Build

Build the Docker image:

```sh
docker build -t gtc_to_vcf .
```

### Usage

Set the full path to the location of the GTC files (`gtc_dir`), manifest csv or bpm (`manifest`), and reference fasta (`ref`). The reference fasta must be indexed with `samtools faidx`, with an accompanying `.fasta.fai` file. 

```sh
gtc_dir=/path/to/gtcs
manifest=/data/projects/cag/iscan/gsa-manifest-clusterfile/GSA-24v3-0_A1.csv
ref=/data/projects/cag/reference-data/reference-genomes/human_g1k_v37.fasta

docker run --rm \
-u $(id -u):$(id -g) \
-v $gtc_dir:/data \
-v $manifest:/tmp/$(basename ${manifest}) \
-v $ref:/tmp/ref.fasta \
-v ${ref}.fai:/tmp/ref.fasta.fai \
gtc_to_vcf --gtc-paths . --manifest-file /tmp/$(basename ${manifest}) --genome-fasta-file /tmp/ref.fasta --output-vcf-path .

```

## License

Copyright 2018 Illumina

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
