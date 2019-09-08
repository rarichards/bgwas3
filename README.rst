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

1. Create a project directory
1. Put all contigs into a directory named "contigs"
2. Create a tsv file with an 'id' column and one or more 'pheno_xxx' columns
3. Run::

   bgwas3 make full

4. View results with::

   open results/plots/index.html

More documentation coming soon.
