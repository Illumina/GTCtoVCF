#!/usr/bin/env bash

if [[ "$#" -lt 1 ]]
then
    echo "usage: bash download_reference.sh output.fa [build_id]"
    exit 1
fi

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    output_file=`readlink -f ${1}`
elif [[ "$OSTYPE" == "darwin"* ]]; then
    output_file=`greadlink -f ${1}`
else
    echo "Error: Unsupported OS ${OSTYPE}"
	exit 1
fi

output_dir=`dirname ${output_file}`
genome_build=${2:-"37"}

# Mitochondrial refseq moved to a different folder, hence the "mt_remote"
if [ ${genome_build} = "37" ]; then
	remote="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
    mt_remote="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/"
elif [ ${genome_build} = "38" ]; then
    remote="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
    mt_remote="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/"
else
	echo "Error: Unsupported genome build ${genome_build}, valid values are 37,38"
	exit 1
fi

temp_dir=`mktemp -d 2>/dev/null || mktemp -d -t ${output_dir}`
pushd ${temp_dir}

for chrom in `seq 1 22` X Y
do
    wget ${remote}/chr${chrom}.fna.gz
done
wget ${mt_remote}/chrMT.fna.gz

build_fa(){
    echo ">${1}" >> ${2}
    gunzip -c chr${1}.fna.gz | grep -v ">" >> "${2}"
}
for chrom in `seq 1 22` X Y
do
    build_fa $chrom $output_file
done
build_fa "MT" $output_file

if hash samtools 2>/dev/null; then
    samtools faidx ${output_file}
else
    echo "Warning: unable to locate samtools, reference file was not indexed"
fi

popd
rm -r ${temp_dir}
