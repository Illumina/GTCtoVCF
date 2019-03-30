#!/usr/bin/env bash

if [[ "$#" -lt 1 ]]
then
    echo "usage: bash download_reference.sh output.fa [build_id]"
    exit 1
fi

output_file=`readlink -f ${1}`
output_dir=`dirname ${output_file}`
genome_build=${2:-"37"}

if [ ${genome_build} = "37" ]; then
	remote="ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq"
elif [ ${genome_build} = "38" ]; then
    remote=ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/
else
	echo "Error: Unsupported genome build ${genome_build}, valid values are 37,38"
	exit 1
fi

temp_dir=`mktemp -d 2>/dev/null || mktemp -d -t ${output_dir}`
pushd ${temp_dir}

for chrom in `seq 1 22` X Y MT
do
    wget ${remote}/*_ref_*chr${chrom}.fa.gz
done

for chrom in `seq 1 22` X Y MT
do
    echo ">${chrom}" >> ${output_file}
    gunzip -c *_ref_*chr${chrom}.fa.gz | grep -v ">" >> "${output_file}"
done

if hash samtools 2>/dev/null; then
    samtools faidx ${output_file}
else
    echo "Warning: unable to locate samtools, reference file was not indexed"
fi

popd
rm -r ${temp_dir}
