import os
from setuptools import setup

HERE = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(HERE, 'README.rst'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='bgwas3',
    version='0.0.1',
    long_description=LONG_DESCRIPTION,
    url='http://github.com/g-r-eg/bgwas3',
    author='Gregory Leeman',
    author_email='g-r-eg@outlook.com',
    license='GNU General Public License v3.0',
    install_requires=[
        "ruffus",
        "cgacore",
        "argparse",
        "docutils"
        ],
    packages=['bgwas3'],
    entry_points={"console_scripts": ['bgwas3=bgwas3.bgwas3:main']},
    include_package_data=True
    )
