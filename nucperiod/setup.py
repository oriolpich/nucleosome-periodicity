from os import path
from setuptools import setup, find_packages


VERSION = "0.1"
DESCRIPTION = "Utilities for Nucleosome Periodicity"

directory = path.dirname(path.abspath(__file__))


# Get the long description from the README file
with open(path.join(directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='nucperiod',
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="https://bitbucket.org/bbglab/nucleosome-periodicity",
    author="Barcelona Biomedical Genomics Lab",
    author_email="bbglab@irbbarcelona.org",
    packages=find_packages(),
    install_requires=['lmfit', 'matplotlib','numpy', 'rpy2', 'pandas'],
)