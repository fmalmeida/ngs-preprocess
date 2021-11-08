.. _config:

Configuration File
******************

To download a configuration file template users just need to run: ``nextflow run fmalmeida/ngs-preprocess [--get_config]``

Using a config file your code is a lot more cleaner and concise: ``nextflow run fmalmeida/ngs-preprocess -c [path-to-config]``

Default configuration:
""""""""""""""""""""""

.. literalinclude:: ../nextflow.config
   :language: groovy