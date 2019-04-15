FROM debian:stretch-slim
WORKDIR /opt/gtc_to_vcf

RUN apt-get update && apt-get install -y wget bzip2
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN printf '\nyes\n' | bash Miniconda2-latest-Linux-x86_64.sh
RUN /root/miniconda2/bin/conda install -c miniconda numpy=1.11.2
RUN /root/miniconda2/bin/conda install -c bioconda pyvcf=0.6.8
RUN /root/miniconda2/bin/conda install -c bioconda pysam=0.9.0

COPY . /opt/gtc_to_vcf

ENTRYPOINT [ "/root/miniconda2/bin/python", "/opt/gtc_to_vcf/gtc_to_vcf.py" ]
