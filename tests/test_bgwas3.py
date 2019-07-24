import pytest
import os
import subprocess
import shutil
import filecmp

test_dir = os.getcwd()
ref_dir = test_dir + "/ref"
script = test_dir + "/../bgwas3/bgwas3.py"
if not os.path.exists("temp"):
    os.mkdir("temp")

def runStep(step_name, req_files, outfile, yml, local=False): # {{{

    os.chdir(test_dir)
    temp_dir = test_dir + "/temp/" + step_name

    if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)

    os.mkdir(temp_dir)

    for f in req_files:
        source = ref_dir + "/" + f
        target = temp_dir + "/" + f
        assert os.path.exists(source), "Reference " +  f + " not found in reference folder."
        if os.path.isfile(source):
            print("Copying reference file '" + f + "' to temp directory:", end=" ")
            shutil.copy(source, target)
            print("done")
        if(os.path.isdir(source)):
            print("Copying reference directory '" + f + "' to temp directory:", end=" ")
            shutil.copytree(source, target)
            print("done")

    with open(temp_dir + "/pipeline.yml", "w") as f:
        f.write(yml)

    os.chdir(temp_dir)

    statement = "python " + script + " make " + step_name

    if local:
        statement += " --local"

    with open("run.sh", "w") as f:
        f.write(statement)

    print("Running: " + statement)
    os.system(statement + " 2>&1 | tee test.log")

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

    out = temp_dir + "/" + outfile
    ref = ref_dir + "/" + outfile
    assert os.path.exists(out), "Output '" + outfile + "' does not exist"
    
    if os.path.isdir(out):
        assert os.listdir(out), "Output directory '" + outfile + "' is empty"
        assert os.listdir(out) == os.listdir(ref), "Output directory '" + outfile + "' does not match reference directory" 

    if os.path.isfile(out):
        assert filecmp.cmp(out, ref), "Output file '" + outfile + "' does not match reference"

# }}}
    
@pytest.mark.cluster
def test_fsm():
    req_files = ["fastq.dir", "contigs.dir"]
    runStep("fsm", req_files, "kmers.gz", "")

@pytest.mark.local
def test_splitPhenos():
    req_files = ["phenos.tsv"]
    runStep("splitPhenos", req_files, "phenos.dir", "", local=True)

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
