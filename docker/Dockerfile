FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    build-essential git wget \
    xsltproc \
    libssl-dev

WORKDIR /
RUN git clone https://github.com/KeccakTeam/KeccakCodePackage
WORKDIR KeccakCodePackage
RUN make generic64/libXKCP.a
RUN cp -a bin/generic64/libXKCP.a /usr/local/lib
RUN cp -a bin/generic64/libXKCP.a.headers /usr/local/include


WORKDIR /
RUN wget -c https://classic.mceliece.org/nist/mceliece-20221023.tar.gz
RUN tar xzf mceliece-20221023.tar.gz

WORKDIR /mceliece-20221023/Reference_Implementation/kem/mceliece348864
RUN sed -i 's/keccak/XKCP/g' crypto_hash.h build
COPY leak_kem.c nist/leak_kem.c
COPY build.leak pk_gen.patch ./
RUN patch < pk_gen.patch
RUN sh build.leak
RUN cp leak /leak-mceliece348864

WORKDIR /mceliece-20221023/Reference_Implementation/kem/mceliece460896
RUN sed -i 's/keccak/XKCP/g' crypto_hash.h build
COPY leak_kem.c nist/leak_kem.c
COPY build.leak pk_gen.patch ./
RUN patch < pk_gen.patch
RUN sh build.leak
RUN cp leak /leak-mceliece460896

WORKDIR /mceliece-20221023/Reference_Implementation/kem/mceliece6688128
RUN sed -i 's/keccak/XKCP/g' crypto_hash.h build
COPY leak_kem.c nist/leak_kem.c
COPY build.leak pk_gen.patch ./
RUN patch < pk_gen.patch
RUN sh build.leak
RUN cp leak /leak-mceliece6688128

WORKDIR /mceliece-20221023/Reference_Implementation/kem/mceliece6960119
RUN sed -i 's/keccak/XKCP/g' crypto_hash.h build
COPY leak_kem.c nist/leak_kem.c
COPY build.leak pk_gen.patch ./
RUN patch < pk_gen.patch
RUN sh build.leak
RUN cp leak /leak-mceliece6960119

WORKDIR /mceliece-20221023/Reference_Implementation/kem/mceliece8192128
RUN sed -i 's/keccak/XKCP/g' crypto_hash.h build
COPY leak_kem.c nist/leak_kem.c
COPY build.leak pk_gen.patch ./
RUN patch < pk_gen.patch
RUN sh build.leak
RUN cp leak /leak-mceliece8192128

WORKDIR /

RUN apt-get update && apt-get install -y \
    python3
RUN wget -c https://botan.randombit.net/releases/Botan-3.1.1.tar.xz
RUN tar xJf Botan-3.1.1.tar.xz
WORKDIR Botan-3.1.1
RUN ./configure.py
COPY botan-pkgen.patch ./
RUN patch -p1 < botan-pkgen.patch
RUN make -j 8
COPY botan /leak-botan

WORKDIR /

