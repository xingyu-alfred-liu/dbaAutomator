import os
from setuptools import setup, find_packages

os.chdir('dbaAutomator/bader')
os.system('make bader')
os.chdir('../../')

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
        name='dbaAutomator',
        version='0.1.0',
        description='Automate Double-Bader Analysis (dba) Process',
        long_description=long_description,
        long_description_content_type="text/markdown",
        author='Xingyu (Alfred) Liu',
        author_email='xingyul1@andrew.cmu.edu, xingyu.alfred.liu@gmail.com',
        url="https://github.com/BLABABA/dbaAutomator.git",
        packages=find_packages(),
        install_requires=['numpy', 'pymatgen', 'ase'],
)