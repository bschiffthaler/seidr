FROM fedora
RUN dnf update -y && dnf install -y gcc gcc-c++ gcc-gfortran cmake git \
    nano tclap zlib-devel wget bzip2-devel bzip2 make \
    csh xz-devel findutils hostname which libicu-devel armadillo-devel gsl-devel \
    glpk-devel openmpi-devel boost-devel

RUN ln -s /usr/lib64/openmpi/bin/* /usr/local/bin/. && \
    ln -s /usr/lib64/openmpi/lib/* /usr/local/lib/.

WORKDIR /build

# Install HTSLIB
WORKDIR /build
RUN wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 && \
    tar -xf htslib-1.8.tar.bz2 && \
    cd htslib-1.8 && \
    ./configure --prefix=/usr && \
    make -j4 && \
    make install 

WORKDIR /build
RUN git clone --recursive https://github.com/bschiffthaler/seidr

WORKDIR /build/seidr

RUN mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make && \
    cp anoverence mi correlation el-ensemble genie3 llr-ensemble narromi pcor plsnet seidr svm-ensemble tigress /usr/local/bin

WORKDIR /
RUN rm -rf /build
