FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing all software requirements for the fmalmeida/ngs-preprocess pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda install -y -c conda-forge mamba
RUN mamba env create --quiet -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/ngs-preprocess-2.3/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ngs-preprocess-2.3 > ngs-preprocess-2.3.yml