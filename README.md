# Steering Patchy Particles using Multivalent Electrolytes

This repository contain a Jupyter Notebook for studying charged, patchy particles
in electrolyte solution using Metropolis Monte Carlo simulations. The layout is as follows:

- `montecarlo.ipynb` - Jupyter Notebook for compiling, running, and plotting all data
- `mc/` - Dicrectory with Monte Carlo data and C++ code
- `environment.yml` - Environment file for setting up dependencies for Anaconda

To open this Notebook, use [Conda](https://www.continuum.io/downloads) and make sure all required packages are loaded
by issuing the following in a terminal:

    conda env create -f environment.yml
    source activate cppm
    jupyter-notebook montecarlo.cpp
