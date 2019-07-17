import pytest
import os
import tempfile

def test_0():
    x = 4
    y = 4
    assert x == y, str(x) + " does not equal " + str(y)

def test_1():
    os.system("python ../bgwas3/bgwas3.py show full > 1.log")
