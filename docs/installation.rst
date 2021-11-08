.. _installation:

Installation
************

Dependencies
============

This pipeline `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ requires `Docker <https://www.docker.com/>`_, `Singularity <https://sylabs.io/singularity/>`_ or `conda <https://conda.io/>`_ to run. Please read the information about how to proper :ref:`choose between conda, docker and singularity profiles<profiles>`.

1. Installing Nextflow

   .. code-block:: bash

      curl -s https://get.nextflow.io | bash

2. Download the pipeline

   .. code-block:: bash

      nextflow pull fmalmeida/ngs-preprocess

3. Test your installation

   .. code-block:: bash

      nextflow run fmalmeida/ngs-preprocess --help

4. Download required tools

   .. code-block:: bash
      
      # for docker
      docker pull fmalmeida/ngs-preprocess:v2.3

      # for singularity
      singularity pull docker://fmalmeida/ngs-preprocess:v2.3

      # for conda
      wget https://github.com/fmalmeida/ngs-preprocess/raw/master/environment.yml
      [mamba|conda] env create -f environment.yml

5. (Optional) Install nf-core

   ``pip install nf-core>=1.10``

.. note::

  Now, everything is set up and ready to run. Remember to always keep your Docker images up to date (Docker pull will always download the latest).
