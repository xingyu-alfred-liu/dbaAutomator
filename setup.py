from setuptools import setup

setup(
        name='dbaAutomator',
        version='1.0',
        description='Automate Double-Bader Analysis (dba) Process',
        author='Xingyu (Alfred) Liu',
        author_email='xingyul1@andrew.cmu.edu, xingyu.alfred.liu@gmail.com',
        packages=['dbaAutomator'],  #same as name
        install_requires=['numpy', 'pymatgen>=2018.11.6', 'ase>=3.17'], #external packages as dependencies
        scripts=[
            'src/functions.py',
            'src/structio.py',
            'src/ref.py',
            'src/core.py',
            ]
)