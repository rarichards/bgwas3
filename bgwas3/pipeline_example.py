import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

from ruffus import *
import os
import sys

# load options from the config file
#P.get_parameters(
#        ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
#            "../pipeline.yml",
#            "pipeline.yml"]
#         )


def echo_this():
    '''print a message'''

    statement = "echo 'This is working'"

    P.run(statement)


@follows(echo_this)
def touch_file():
    ''' dummy function touch file'''

    statement = 'touch test_file.txt'

    P.run(statement)

@follows(touch_file)
def full():
    '''run all'''
    pass

# Finish and exit:
if __name__ == '__main__':
    sys.exit(P.main(sys.argv))
