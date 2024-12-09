From ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive

# RUN apt-get --allow-releaseinfo-change update
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils 
RUN apt-get install -f --fix-broken -y --no-install-recommends \
  locales \
  python3-pip \
  ca-certificates \
  netbase \
  gnupg \
  dirmngr \
  git \
  openssh-client \
  rpm \
  sshpass \
  valgrind \
  gdb \
  clang-format \
  ninja-build \ 
  doxygen \
  graphviz \
  curl \
  wget \
  build-essential \
  cmake \
  libboost-all-dev \
  pkg-config \
  make \
  cmake \
  g++ \
  libeigen3-dev \
  gfortran \
  libmetis-dev \
  libopenvdb-dev \
  libboost-all-dev \
  libtbb2 \
  libtbb-dev \
  libz-dev \
  libgmp3-dev \
  libmpc-dev \
  vim \
  libeigen3-dev \
  libcgal-dev \
  libblas-dev liblapack-dev libsuitesparse-dev \
  xorg-dev libglu1-mesa-dev \
  libglu1-mesa-dev freeglut3-dev mesa-common-dev \
  libgmsh-dev gmsh\
  libfmt-dev libglm-dev libopenblas-dev liblapack-dev liblapack3 libopenblas-base \
  && apt-get clean

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  curl ca-certificates gpg-agent software-properties-common && \
  rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install mesa-utils -y --no-install-recommends  && rm -rf /var/lib/apt/lists/*
# repository to install Intel(R) oneAPI Libraries
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | gpg --dearmor | tee /usr/share/keyrings/intel-oneapi-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/intel-oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main " > /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  curl ca-certificates gpg-agent software-properties-common && \
  rm -rf /var/lib/apt/lists/*
# repository to install Intel(R) GPU drivers
RUN curl -fsSL https://repositories.intel.com/graphics/intel-graphics.key | gpg --dearmor | tee /usr/share/keyrings/intel-graphics-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/intel-graphics-archive-keyring.gpg arch=amd64] https://repositories.intel.com/graphics/ubuntu jammy flex" > /etc/apt/sources.list.d/intel-graphics.list

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  ca-certificates build-essential pkg-config gnupg libarchive13 openssh-server openssh-client wget net-tools git intel-oneapi-mkl-devel-2024.2 intel-oneapi-tbb-devel-2022.0 && \
  rm -rf /var/lib/apt/lists/*

RUN cp -r /usr/include/eigen3/* /usr/local/include/

RUN git clone --recurse-submodules https://github.com/lnumerow-ethz/geogram && cd geogram && mkdir build && cd build && cmake .. && make -j8 install && cd ../..

ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES all
