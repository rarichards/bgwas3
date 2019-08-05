from setuptools import setup
import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

install_requires = [
    "ruffus",
    "cgacore",
    "argparse"
    ]

setup(
    name = 'bgwas3',
    version = '0.0.1',
    long_description = long_description,
    url = 'http://github.com/g-r-eg/bgwas3',
    author = 'Gregory Leeman',
    author_email = 'g-r-eg@outlook.com',
    license = 'GNU General Public License v3.0',
    packages = ['bgwas3'],
    entry_points = {"console_scripts": ['bgwas3 = bgwas3.bgwas3:main']},
    install_requires = [
        "ruffus",
        "cgacore",
        "argparse"
        ]
    )
