from setuptools import setup

setup(
    name='bgwas3',
    version='1.0',
    url='http://github.com/g-r-eg/bgwas3',
    author='Gregory Leeman',
    author_email='g-r-eg@outlook.com',
    license='GNU General Public License v3.0',
    packages=['bgwas3'],
    entry_points={"console_scripts": ['bgwas3=bgwas3.bgwas3:main']}
)
