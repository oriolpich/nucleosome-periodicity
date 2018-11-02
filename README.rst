
README
======

This folder contains the code used for `10.1016/j.cell.2018.10.004 <https://doi.org/10.1016/j.cell.2018.10.004>`_.
If you use this code in a publication, please cite:

.. admonition:: Citation
   :class: note

   Oriol Pich, Ferran Muiños, Radhakrishnan Sabarinathan, Iker Reyes-Salazar, Abel Gonzalez-Perez,
   Nuria Lopez-Bigas, Somatic and germline mutation periodicity follow the orientation of the DNA minor
   groove around nucleosomes, Cell (2018) doi: `10.1016/j.cell.2018.10.004 <https://doi.org/10.1016/j.cell.2018.10.004>`_



A brief description of the structure of this repo:

- `accessibility <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/accessibility/accessibility.ipynb>`_: accessibility data analysis
- `ancestral <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/ancestral/ancestral.ipynb>`_: ancestral states analysis
- `damage <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/damage/damage.ipynb>`_: damage data (UV and NMP) analysis
- `figures <https://bitbucket.org/bbglab/nucleosome-periodicity/src/master/figures/>`_: code to generate the figures and tables for the paper
- `germline <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/germline/germline.ipynb>`_: germline data analysis
- `increase <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/increase/increase.ipynb>`_: code for the increase of mutation rate analysis
- `mutations <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/mutations/mutations.ipynb>`_:  mutational data analysis
- `nucleosomes <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/nucleosomes/nucleosomes.ipynb>`_:  computation of the dyads positions
- `periodicity <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/periodicity/periodicity.ipynb>`_:  WW periodicity analysis
- `rotational <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/rotational/rotational.ipynb>`_:  rotational classification of the nucleosomes
- `signatures <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/signatures/signatures.ipynb>`_:  analysis of the signatures of the mutational data
- `simulation <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/simulation/simulation.ipynb>`_:  simulation

Each folder contains a notebook with a brief description and
the requirements (notebooks that need to be executed).


Running this software
---------------------

These analysis have been perform using software in Python, R and GNU bash.

We have created a set of `Jupyter notebooks <http://jupyter.org/>`_
that you can run if you are interested in re-running partially or
totally our analysis.
In each notebook you will find further details for running them.

Requirements
************

To be able to run those notebooks you need to have the following
software installed (we also indicate the version so you can
reproduce the exact same results):

   Python (3.5.6) Packages:

   - `ipykernel <https://pypi.org/project/ipykernel/>`_ (4.8.2)
   - `numpy <http://www.numpy.org/>`_ (1.15.1)
   - `pandas <https://pandas.pydata.org/>`_ (0.23.4)
   - `scipy <https://www.scipy.org/>`_ (1.1.0)
   - `matplotlib <https://matplotlib.org/>`_ (2.2.3)
   - `statsmodels <https://www.statsmodels.org/stable/index.html>`_ (0.9.0)
   - `click <http://click.pocoo.org>`_ (6.7)
   - `tqdm <https://pypi.org/project/tqdm>`_ (4.25.0)
   - `intervaltree <https://pypi.org/project/intervaltree>`_ (2.1.0)
   - `lmfit <https://lmfit.github.io/lmfit-py>`_ (0.9.11)
   - `rpy2 <https://rpy2.readthedocs.io/en/latest/>`_ (2.7.8)
   - `bgreference <https://bitbucket.org/bgframework/bgreference>`_ (0.5)
   - `xlrd <http://www.python-excel.org/>`_ (1.1.0)

   Python (2.7.15) Packages:

   - `crossmap <http://crossmap.sourceforge.net>`_ (0.2.7) [#envcrossmap]_

   R (3.4.3) packages:

   - `deconstructSigs <https://github.com/raerose01/deconstructSigs>`_ (1.8.0) [#envdeconstruct]_
   - `sigfit <https://github.com/kgori/sigfit>`_ (1.0.0)  [#envsigfit]_

   Other software:

   - `bedops <https://bedops.readthedocs.io/en/latest/>`_ (2.4.35)
   - `bedgraphtobigwig <http://hgdownload.soe.ucsc.edu/admin/exe/>`_ (366)
   - `bigwigtowig <http://hgdownload.soe.ucsc.edu/admin/exe/>`_ (366)
   - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_ (v2.27.1)
   - `sra-tools <https://github.com/ncbi/sra-tools>`_ (2.9.1)
   - `bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_ (1.2.2)
   - `vcftools <https://vcftools.github.io/index.html>`_ (0.1.16)
   - `fa-split <http://hgdownload.soe.ucsc.edu/admin/exe/>`_ (366)
   - `bwtool <https://github.com/CRG-Barcelona/bwtool/wiki>`_ (1.0) [#noconda]_
   - `awk <http://www.cs.princeton.edu/~bwk/btl.mirror/>`_ (4.0.2) [#noconda]_
   - `sort <http://www.gnu.org/software/coreutils/>`_ (8.22) [#noconda]_
   - `gzip <https://www.gnu.org/software/gzip/>`_ (1.5) [#noconda]_
   - `grep <https://www.gnu.org/software/grep/manual/grep.html>`_ (2.20) [#noconda]_
   - `basename <http://www.gnu.org/software/coreutils/>`_ (8.22) [#noconda]_
   - `mktemp <http://www.gnu.org/software/coreutils/>`_ (8.22) [#noconda]_

In addition, we have created a Python package named ``nucperiod`` that contains a set of
python scripts that we have used during our analysis.
In can be installed with pip:

.. code:: python

   cd nucperiod
   pip install .

For some of the analyses (those where CrossMap, DeconstructSigs and SigFit are involved)
we already prepared three `conda environments <https://conda.io/docs/>`_:

- ``env_crossmap`` environment for Crossmap as it is a Python 2.7 tool
  (you can create it with ``env_crossmap.yml``)
- ``env_deconstructsigs`` environment for the deconstructSigs
  R package (use the ``env_deconstructsigs.yml`` to replicate it)
- ``env_nucperiod_sigfit`` environment for the SigFit R package.
  Please, note that the package is not installed in that environment
  and you need to install it manually


Notes
*****

Most of the intermediate files generated while running any notebook
are most likely not used for further analysis.
However, we have decided not to remove them so you can
check them if needed.

Compressing most of the files is not needed, however, we
have decided to do that in order to save disk space.

The scripts that you can find in the ``scripts`` directories
are documented for further info.
If you want to check which parameters
each script accepts, use the ``--help`` flag
(`python <script> --help`).


Fixing datasets versions
------------------------

This project makes use of datasets available thought the
`bgdata <https://bitbucket.org/bgframework/bgdata>`_.
This package will try to download the latest version,
however, you can fix the version of these datasets easily.
After installing the package, update the file
``~/.bbglab/bgdata.conf`` to add the following lines::

    [datasets/genomereference/hg19]
    build = 20150724
    [datasets/genomereference/tair10]
    build = 20180810
    [datasets/genomereference/saccer3]
    build = 20180720
    [datasets/genomereference/dm3]
    build = 20180904
    [datasets/genomereference/mm9]
    build = 20171103

----

.. [#noconda] This software was *not* installed within
   a conda environment.

.. [#envcrossmap] This package has been installed in a separate environment
   named as ``env_crossmap``

.. [#envdeconstruct] This package has been installed in a separate environment
   named as ``env_deconstructsigs``

.. [#envsigfit] This package has been installed in a separate environment
   named as ``env_sigfit``

