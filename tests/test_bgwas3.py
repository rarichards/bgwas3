import pytest
import os
import subprocess
import shutil
import filecmp

def test_listContigs():
    if os.path.exists("temp") and os.path.isdir("temp"):
        shutil.rmtree("temp")

    shutil.copytree("ref/contigs", "temp/contigs")
    shutil.copytree("ref/fastq", "temp/fastq")
    shutil.copy("ref/phenos.tsv", "temp/phenos.tsv")

    shutil.copy("ref/pipeline.yml", "temp/pipeline.yml")

    os.chdir("temp")
    statement = "python ../../bgwas3/bgwas3.py make listContigs --local"
    os.system(statement)
    os.chdir("../")

    assert os.path.isfile("temp/contig_list.txt"), "Ouput file does not exist"
    assert filecmp.cmp("temp/contig_list.txt", 'ref/contig_list.txt'), "Output file does not match reference"

    shutil.rmtree("temp")

def test_getKmers():
    if os.path.exists("temp") and os.path.isdir("temp"):
        shutil.rmtree("temp")

    shutil.copytree("ref/contigs", "temp/contigs")
    shutil.copytree("ref/fastq", "temp/fastq")
    shutil.copy("ref/phenos.tsv", "temp/phenos.tsv")

    shutil.copy("ref/pipeline.yml", "temp/pipeline.yml")

    os.chdir("temp")
    statement = "python ../../bgwas3/bgwas3.py make getKmers --local"
    os.system(statement)
    os.chdir("../")

    assert os.path.isfile("temp/kmers.txt.gz"), "Ouput file does not exist"

    shutil.rmtree("temp")
