from setuptools import setup
from codecs import open
from os import path
import os
import re
import io

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

install_requires = []

with open(os.path.join(here, 'requirements.txt'), encoding='utf-8') as required:
    for line in (required):
        if not line.startswith('#') and not line.startswith('\n'):
            line = line.strip()
            install_requires.append(line)

setup(
    name = 'bgwas3',
    version = '0.1',
    description = 'Sequence Elements Enrichment Analysis (SEER), python implementation',
    long_description = long_description,
    url = 'http://github.com/g-r-eg/bgwas3',
    author = 'Gregory Leeman',
    author_email = 'gregoryleeman@outlook.com',
    license = 'GNU General Public License v3.0',
    packages = ['bgwas3'],
    entry_points = {
        "console_scripts": [
            'bgwas3 = bgwas3.__main__:main',
            ]
    },
    install_requires = install_requires,
    zip_safe = False
    )

