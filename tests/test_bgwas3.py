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

def runStep(step_name, out, local=False): # {{{

    # make temporary directory with a copy of all reference files {{{
    os.chdir(test_dir)
    temp_dir = test_dir + "/temp/" + step_name

    if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)

    shutil.copytree(ref_dir, temp_dir)

    os.chdir(temp_dir)

    # }}}

    # make all pipeline tasks up to date with touch {{{
    statement = "python " + script + " touch " + step_name
    assert subprocess.check_call(statement, shell=True) == 0, "Can't touch pipeline"

    # }}}

    # remove output to be regenerated {{{
    for f in out:
        assert os.path.exists(f), "Reference " +  f + " not found in reference folder."
        if os.path.isdir(f):
            shutil.rmtree(f)
        elif os.path.isfile(f):
            os.remove(f)

    print(os.listdir())

    # }}}

    # make statment (force up to date step of interest)
    statement = "python " + script + " make " + step_name + " -f " + step_name

    # run statement
    if local:
        statement += " --local"
        assert subprocess.check_call(statement, shell=True) == 0, "Can't run step"
    else:
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

# }}}

# }}}
    
@pytest.mark.cluster
def test_fsm():
    runStep("fsm", ["kmers.gz"])

@pytest.mark.local
def test_splitPhenos():
    runStep("splitPhenos", ["phenos"], local=True)

@pytest.mark.local
def test_prokka():
    runStep("prokka", ["prokka"], local=True)

@pytest.mark.cluster
def test_roary():
    runStep("roary")

@pytest.mark.local
def test_distanceFromTree():
    runStep("distanceFromTree", ["distances.tsv"], local=True)

@pytest.mark.cluster
def test_pyseer():
    runStep("pyseer", ["pyseer"])

@pytest.mark.local
def test_pyseer_local():
    runStep("pyseer", ["pyseer"], local=True)
<<<<<<< HEAD

@pytest.mark.local
def test_makeRefList():
    runStep("makeRefList", ["ref.txt"], local=True)

@pytest.mark.cluster
def test_mapKmers():
    runStep("mapKmers", ["maps"])
=======
>>>>>>> 364210e25067eac66f5cc0e10e1400969c251f5a
