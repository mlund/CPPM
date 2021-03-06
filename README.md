# Steering Patchy Particles using Multivalent Electrolytes

_Alexei Abrikossov, Björn Stenqvist, Mikael Lund, 2017._

[DOI:10.1039/C7SM00470B](http://dx.doi.org/10.1039/C7SM00470B).

This repository contain a [Jupyter](http://jupyter.org) Notebook for studying charged, patchy particles
in electrolyte solution using Metropolis Monte Carlo simulations. The layout is as follows:

- `montecarlo.ipynb` - Jupyter Notebook for compiling, running, and plotting all simulation data
- `environment.yml` - Environment file for setting up dependencies for Anaconda
- `mc/` - Directory with Monte Carlo data and C++ code
- `md/` - Directory with Notebook for MD simulations using OpenMM (experimental!)

To open this Notebook, install python via [(Mini)conda](https://www.continuum.io/downloads) and make sure all required packages are loaded
by issuing the following terminal commands,

    conda env create -f environment.yml
    source activate cppm
    jupyter-notebook montecarlo.ipynb
