import pytest
import os
import tempfile

ref_dir = os.path.abspath(os.path.join('tests', 'dat'))
test_dir = tempfile.mkdtemp()
os.chdir(test_dir)

def test_0():
    x = 4
    y = 4
    assert x == y, str(x) + " does not equal " + str(y)
