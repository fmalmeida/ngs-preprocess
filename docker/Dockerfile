FROM ubuntu:18.04

MAINTAINER Felipe Marques <180039822@aluno.unb.br>
SHELL ["/bin/bash", "-c"]

RUN apt-get update

RUN apt-get upgrade -y

RUN apt-get install -y python python3 python-pip git default-jre fastqc wget cutadapt \
	curl python-pkg-resources python3-pip python3-pkg-resources git python-setuptools \
	python-dev cython libhdf5-serial-dev

# INSTALL MINICONDA
RUN apt-get update && apt-get upgrade -y && apt-get install -y git libc6-dev build-essential wget && \
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
	bash /miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# Update do latest
RUN conda update -y -n base conda

# Configure Miniconda
RUN conda config --add channels defaults && \
	conda config --add channels bioconda && \
	conda config --add channels conda-forge && \
	conda init --all

RUN pip3 install --user --upgrade cutadapt

RUN pip install --upgrade pip ; pip3 install --upgrade pip

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y python3-matplotlib python-matplotlib

# INSTALL PORETOOLS
#RUN git clone https://github.com/arq5x/poretools ; cd poretools ; python2.7 setup.py install
#RUN conda create -y -n poretools && conda install -y -n poretools -c bioconda poretools

# Install Trim Galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.tar.gz -o trim_galore.tar.gz ; tar xvzf trim_galore.tar.gz ; ln -rs TrimGalore-0.6.0/* /usr/local/bin

# INSTALL PORECHOP
RUN git clone https://github.com/rrwick/Porechop.git ; cd Porechop ; python3 setup.py install

# INSTALL LIGHTER ILLUMINA READ CORRECTION
RUN git clone https://github.com/mourisl/Lighter.git ; cd Lighter/ ; make ; ln -rs lighter /usr/local/bin

# INSTALL PEAR
#RUN conda install -y pear ; conda update -y pear ; conda install -y h5py

# INSTALL FLASH
RUN conda create -y -n flash ; conda install -y -n flash -c bioconda flash

# INSTALL LoRDEC hybrid long reads correction
#RUN wget http://gatb-tools.gforge.inria.fr/versions/bin/LoRDEC-0.5.3-Linux.tar.gz ; tar zxvf LoRDEC-0.5.3-Linux.tar.gz ; ln -rs LoRDEC-0.5.3-Linux/bin/* /usr/local/bin

# INSTALL lorma longreads only reads correction
#RUN conda install -c atgc-montpellier lorma

# INSTALL NANOPOLISH
RUN apt-get install -y libhdf5-dev ; conda config --add channels defaults ; conda config --add channels conda-forge ; conda config --add channels bioconda ; conda install -y nanopolish pear

# INSTALL TOOLS FOR PACBIO READS EXTRACTION
RUN apt-get install -y pbh5tools
RUN conda create -y -n pbtools && conda install -y -n pbtools -c bioconda bax2bam bam2fastx lima

# INSTALL FMLRC
#RUN git clone https://github.com/holtjma/fmlrc.git ; cd fmlrc ; make ; ln -rs fmlrc /usr/local/bin ; ln -rs fmlrc-convert /usr/local/bin

# INSTALL ROPEBWT2
#RUN git clone https://github.com/lh3/ropebwt2.git ; cd ropebwt2 ; make ; ln -rs ropebwt2 /usr/local/bin

# INSTALL NANOPACK
RUN conda create -y -n nanopack && conda install -y -n nanopack -c bioconda nanoplot nanostat nanoqc

# INSTALL FILTLONG
#RUN conda create -y -n filtlong && conda install -y -n filtlong -c bioconda filtlong

# INSTALL pycoQC
RUN conda create -y -n pycoQC python=3.6 && \
  conda install -y -n pycoQC -c anaconda pip && \
  source activate pycoQC && \
  python3 -m pip install pycoQC

WORKDIR /work
