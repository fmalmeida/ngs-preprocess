FROM ubuntu:18.04

LABEL MAINTAINER "Felipe Marques <180039822@aluno.unb.br>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update

RUN apt-get upgrade -y

RUN apt-get install -y python python3 python-pip git default-jre fastqc wget \
			curl python-pkg-resources python3-pip python3-pkg-resources git python-setuptools \
			python-dev cython libhdf5-serial-dev

# INSTALL MINICONDA
RUN apt-get update && apt-get upgrade -y && apt-get install -y git libc6-dev build-essential wget && \
			wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
			bash /miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# Update do latest
RUN conda update -y -n base conda && \
			conda clean -afy

# Configure Miniconda
RUN conda config --add channels defaults && \
			conda config --add channels bioconda && \
			conda config --add channels conda-forge && \
			conda init --all

RUN conda install -y -c bioconda cutadapt && \
			cutadapt --version

RUN pip install --upgrade pip && \
			pip3 install --upgrade pip

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y python3-matplotlib python-matplotlib

# Install Trim Galore
RUN export TRIM_GALORE_TAG=$(curl -s "https://api.github.com/repos/FelixKrueger/TrimGalore/releases/latest" | grep "tag_name" | cut -d ":" -f 2,3  | tr -d '"' | tr -d "," | tr -d " ") && \
			curl -fsSL "https://github.com/FelixKrueger/TrimGalore/archive/${TRIM_GALORE_TAG}.tar.gz" -o trim_galore.tar.gz && \
			tar xvzf trim_galore.tar.gz ; ln -rs TrimGalore-*/* /usr/local/bin

# INSTALL PORECHOP
RUN git clone https://github.com/rrwick/Porechop.git ; cd Porechop ; python3 setup.py install

# INSTALL LIGHTER ILLUMINA READ CORRECTION
# RUN git clone https://github.com/mourisl/Lighter.git ; cd Lighter/ ; make ; ln -rs lighter /usr/local/bin
RUN conda install -c bioconda lighter nomkl && \
			conda clean -afy

# INSTALL FLASH
RUN conda create -y -n flash && \
			conda install -y --freeze-installed -n flash -c bioconda flash nomkl && \
			conda clean -afy

# INSTALL TOOLS FOR PACBIO READS EXTRACTION
RUN apt-get install -y pbh5tools
RUN conda create -y -n pbtools-h5 && \
			conda install -y --freeze-installed -n pbtools-h5 -c bioconda 'bax2bam=0.0.9' 'pbbam=0.18.0' nomkl && \
			conda clean -afy
RUN conda create -y -n pbtools && \
			conda install -y -n pbtools -c bioconda bam2fastx pbccs lima 'pbbam>=0.23.0' nomkl && \
			conda clean -afy			

# INSTALL NANOPACK
RUN conda create -y -n nanopack python=3.7 && \
			source activate nanopack && \
			pip install nanopack biopython==1.78

# INSTALL FILTLONG
RUN conda create -y -n filtlong && \
			conda install -y --freeze-installed -n filtlong -c bioconda filtlong nomkl && \
			conda clean -afy

# INSTALL pycoQC
RUN conda create -y -n pycoQC python=3.6 && \
			source activate pycoQC && \
			conda install -y --freeze-installed -c aleg pycoqc nomkl && \
			conda config --add channels r && \
			conda config --add channels bioconda && \
			conda install -y --freeze-installed pysam nomkl

# Fix matplotlib
RUN source activate nanopack && \
			pip uninstall -y matplotlib && \
			pip install matplotlib==3.1.3

# Clean
RUN conda clean -afy \
			&& find /miniconda/ -follow -type f -name '*.a' -delete \
			&& find /miniconda/ -follow -type f -name '*.pyc' -delete \
			&& find /miniconda/ -follow -type f -name '*.js.map' -delete

WORKDIR /work
