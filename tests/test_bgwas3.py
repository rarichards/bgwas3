import pytest
import os
import subprocess
import shutil
import filecmp
from pathlib import Path

test_dir = os.getcwd()
ref_dir = test_dir + "/ref"
script = test_dir + "/../bgwas3/bgwas3.py"
if not os.path.exists("temp"):
    os.mkdir("temp")

def runStep(step_name, local=False): # {{{

    # make temporary directory 
    os.chdir(test_dir)
    temp_dir = test_dir + "/temp/" + step_name

    if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)

    shutil.copytree(ref_dir, temp_dir)

    #for f in req_files:
    #    source = ref_dir + "/" + f
    #    target = temp_dir + "/" + f
    #    assert os.path.exists(source), "Reference " +  f + " not found in reference folder."
    #    if os.path.isfile(source):
    #        print("Copying reference file '" + f + "' to temp directory:", end=" ")
    #        shutil.copy(source, target)
    #        print("done")
    #    if(os.path.isdir(source)):
    #        print("Copying reference directory '" + f + "' to temp directory:", end=" ")
    #        shutil.copytree(source, target)
    #        print("done")

    #with open(temp_dir + "/pipeline.yml", "w") as f:
    #    f.write(yml)

    os.chdir(temp_dir)

    # make all pipeline tasks up to date
    
    assert subprocess.call("python " + script + " touch full", shell=True) == 0, "Cannot touch the pipeline"
    os.remove("pipeline.log")

    # make statment (force up to date step of interest)
    statement = "python " + script + " make " + step_name + " -f " + step_name

    if local:
        statement += " --local"

    # save statement for debugging later
    with open("statement.sh", "w") as f:
        f.write(statement)

    # run statement
    subprocess.call(statement, shell=True)

    # if cluster fails (dumps core.xxxx) then make a qsub file to test manually instead
    if len([f for f in os.listdir() if f.startswith("core")]) != 0:

        f_name = [f for f in os.listdir() if f.endswith(".sh") and f.startswith("ct")]
        f = open(f_name[0], "r")
        pbs_template = open(test_dir + "/pbs_temp.sh", "r")
        pbs_name = "pbs_" + step_name + ".sh"
        pbs_f = open(pbs_name ,'w')
        lines = pbs_template.readlines() + f.readlines()
        for line in lines:
            pbs_f.write(line)
        f.close()
        pbs_template.close()
        pbs_f.close()
        job_id = os.popen("qsub " + pbs_name).read()
        assert False, "Core dumped (cluster not working). " + pbs_name + " submitted with qsub instead instead (" + job_id + ")"

    # test the contents of the temp directory and the reference directory
    assert filecmp.dircmp(temp_dir, ref_dir), "Output does not match reference"

    #out = temp_dir + "/" + outfile
    #ref = ref_dir + "/" + outfile
    #assert os.path.exists(out), "Output '" + outfile + "' does not exist"
    
    #if os.path.isdir(out):
    #    assert os.listdir(out), "Output directory '" + outfile + "' is empty"
    #    assert os.listdir(out) == os.listdir(ref), "Temporary directory '" + outfile + "' does not match reference" 

    #if os.path.isfile(out):
    #    assert filecmp.cmp(out, ref), "Temporary file '" + outfile + "' does not match reference"

# }}}
    
@pytest.mark.cluster
def test_fsm():
    req_files = ["fastq.dir", "contigs.dir"]
    runStep("fsm", req_files, "kmers.gz", "")

@pytest.mark.local
def test_splitPhenos():
    runStep("splitPhenos", local=True)

@pytest.mark.cluster
def test_prokka():
    req_files = ["fastq.dir", "contigs.dir"]
    runStep("prokka", req_files, "annotations.dir", "", local=True)

@pytest.mark.cluster
def test_roary():
    req_files = ["fastq.dir", "contigs.dir", "annotations.dir"]
    runStep("roary", req_files, "roary.dir", "")

@pytest.mark.local
def test_distanceFromTree():
    req_files = ["phenos.tsv", "phenos.dir"]
    runStep("distanceFromTree", req_files, "distance.tsv", "", local=True)

@pytest.mark.cluster
def test_pyseer():
    req_files = ["fastq.dir", "contigs.dir", "annotations.dir", "roary.dir", "phenos.tsv", "phenos.dir", "distances.tsv"]
    runStep("pyseer", req_files, "annotations.dir", "")
