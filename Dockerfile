# Stage 1: Compile & build wheels for Python deps
FROM alpine:3.11
RUN apk update && apk add --no-cache make git py2-pip python2-dev gcc g++ zlib-dev libcurl curl-dev cython libc-dev ncurses-dev xz-dev bzip2-dev linux-headers
RUN ln -s /usr/include/locale.h /usr/include/xlocale.h
RUN pip install wheel
RUN pip install PyVCF==0.6.8 
RUN pip install pysam==0.11.0 
RUN pip install numpy==1.11.2 
RUN pip wheel PyVCF pysam numpy
RUN mkdir /src && mv *.whl /src

# Step 2: Install Python deps from wheels, add code and tests
FROM alpine:3.11
ADD scripts/* ./*.py /usr/local/bin/
ADD tests /tests
COPY --from=0 /src/*whl /src/
RUN apk update && apk add --no-cache py2-pip python2 xz-dev libcurl
RUN pip install /src/*.whl && rm /src/*whl
RUN mkdir /data && chmod 777 /data

# Run
WORKDIR /data
ENTRYPOINT ["gtc_to_vcf.py"]
