FROM ubuntu:20.04

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    curl \
    build-essential \
    zlib1g-dev 

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN /bin/bash -c "curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-\$(uname -m).sh > mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh"

RUN /bin/bash -c "mamba create -q -y -c bioconda -n gaftools poetry cython==3.0.7 python==3.10 setuptools"
RUN echo "source activate gaftools" > ~/.bashrc
ENV PATH /opt/conda/envs/gaftools/bin:${PATH}
COPY . .
RUN conda init bash && . ~/.bashrc && conda activate gaftools && poetry install
