Quick Start
===========

Full pipeline from config
-------------------------

Create a YAML config file and run the full pipeline:

.. code-block:: bash

   ncountr init > my_analysis.yaml
   # edit my_analysis.yaml with your settings
   ncountr run my_analysis.yaml

Python API
----------

.. code-block:: python

   import ncountr

   # Parse RCC files
   exp = ncountr.read_rcc("path/to/rcc_files/")
   print(exp)  # NanostringExperiment(785 genes x 47 samples, ...)

   # QC
   qc_results = ncountr.qc(exp)

   # Normalize
   ncountr.normalize(exp, method="pos_hk")

   # Differential expression
   de_results = ncountr.de(exp, group_a=["S1", "S2"], group_b=["S3", "S4"])

   # Export to AnnData (for scverse ecosystem)
   adata = ncountr.to_anndata(exp)

Download data from GEO
----------------------

.. code-block:: bash

   ncountr fetch-geo GSE275334 -o data/

Or from Python:

.. code-block:: python

   from ncountr.io.geo import fetch_geo
   rcc_dir = fetch_geo("GSE275334", output_dir="data/")
