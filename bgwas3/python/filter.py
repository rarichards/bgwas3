import argparse

parser = argparse.ArgumentParser(description="Filter associated kmers with bonferoni and generate a stats file")

parser.add_argument('patterns', help='File of patterns from pyseer')
parser.add_argument('assoc', help='File of patterns from pyseer')
parser.add_argument('--alpha', default=0.05, type=float, help='Family-wise error rate')


import sys
import subprocess
from decimal import Decimal

command = "LC_ALL=C sort -u "
if options.cores > 1:
    command +=  "--parallel=" + str(options.cores)
command += (" -S " + str(int(options.memory) - mem_adjust) + "M" +
           " -T " + options.temp +
           " " + options.patterns +
           " | wc -l")

p = subprocess.check_output(command, shell=True, universal_newlines=True)
print("file\tpatterns\tp_thresh")
print(options.patterns + "\t" + p.rstrip() + "\t" + '%.2E' % Decimal(options.alpha/float(p.rstrip())))
