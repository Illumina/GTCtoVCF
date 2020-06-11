#!/usr/bin/env python
import setuptools
from setuptools_scm import get_version


with open("README.md", "r") as fh:
    long_description = fh.read()

    
version = get_version()


setuptools.setup(
    name="gtc_to_vcf",  
#     version=version,
    author="Illumina, Inc.",
    description="Generate a VCF file from Illumina microarray GTC files",
    long_description=long_description,
    #     long_description_content_type="text/markdown",
    url="https://github.com/Illumina/GTCtoVCF",
#     download_url="https://git.illumina.com/Bioinformatics/consensus_sequence/archive/v{}.tar.gz".format(version),
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    scripts=['gtc_to_vcf.py'],
    python_requires='>=3.6',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    requires=['setuptools_scm'],
    install_requires=['setuptools_scm']
)
