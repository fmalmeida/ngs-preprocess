.. _profiles:

Selecting between profiles
**************************

Nextflow profiles are a set of "sensible defaults" for the resource requirements of each of the steps in the workflow, that can be enabled with the command line flag ``-profile``. You can learn more about nextflow profiles at:

+ https://nf-co.re/usage/configuration#basic-configuration-profiles
+ https://www.nextflow.io/docs/latest/config.html#config-profiles

The pipeline have "standard profiles" set to run the workflows with either conda, docker or singularity using the `local executor <https://www.nextflow.io/docs/latest/executor.html>`_, which is nextflow's default and basically runs the pipeline processes in the computer where Nextflow is launched. If you need to run the pipeline using another executor such as sge, lsf, slurm, etc. you can take a look at `nextflow's manual page <https://www.nextflow.io/docs/latest/executor.html>`_ to proper configure one in a new custom profile set in your personal copy of `ngs-preprocess config file <https://github.com/fmalmeida/ngs-preprocess/blob/master/nextflow.config>`_ and take advantage that nextflow allows multiple profiles to be used at once, e.g. ``-profile conda,sge``.

By default, if no profile is chosen, the pipeline will try to load tools from the local machine $PATH. Available pre-set profiles for this pipeline are: ``docker/conda/singularity``, you can choose between them as follows:

* conda

.. code-block:: bash

    # must be executed from the base environment
    # This tells nextflow to load the available ngs-preprocess environment when required
    nextflow run fmalmeida/ngs-preprocess -profile conda [options]

* docker
    
.. code-block:: bash

    nextflow run fmalmeida/ngs-preprocess -profile docker [options]

* singularity
    
.. code-block:: bash

    nextflow run fmalmeida/ngs-preprocess -profile singularity [options]

.. note::
    
    Please use conda as last resource since the packages will not be "frozen and pre-installed", problems may arise, and nextflow will trigger an installation every time which may consume plenty of time.