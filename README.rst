.. image:: https://readthedocs.org/projects/bgwas3/badge/?version=latest
   :target: https://bgwas3.readthedocs.io/en/latest/

.. image:: https://travis-ci.com/g-r-eg/bgwas3.svg?branch=master
   :target: https://travis-ci.com/g-r-eg/bgwas3

bgwas3
======

A pipeline for the pangenome-wide assocation testing of multiple bacteria phenotypes.

Install
-------

::

   conda install bgwas3


Quick Usage
-----------

1. Create and cd into an empty project directory

2. Put all contigs into a sub-directory named "contigs"

3. Create a tsv file with an 'id' column and one or more 'pheno_xxx' columns

4. Run:

::

   bgwas3 make full

5. View results with:

::

   open results/plots/index.html

More documentation coming soon.
