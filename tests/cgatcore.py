import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

from ruffus import *
import os
import sys

def echo_this():
    statement = "echo 'This is working'"
    P.run(statement)

@follows(echo_this)
def touch_file():
    statement = 'touch test_file.txt'
    P.run(statement)

@follows(touch_file)
def full():
    pass

if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
