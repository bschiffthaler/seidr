FROM bschiffthaler/seidr-build

WORKDIR /build

RUN git clone --depth 10 https://github.com/bschiffthaler/seidr

WORKDIR /build/seidr

RUN git submodule sync && git submodule update --init --depth 50

RUN mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j4 && \
    cp anoverence mi correlation el-ensemble genie3 llr-ensemble narromi pcor plsnet seidr svm-ensemble tigress /usr/local/bin

WORKDIR /

RUN rm -rf /build
