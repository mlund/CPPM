# Steering Patchy Particles using Multivalent Electrolytes

This repository contain a Jupyter Notebook for studying charged, patchy particles
in electrolyte solution using Metropolis Monte Carlo simulations. The layout is as follows:

- `montecarlo.ipynb` - Jupyter Notebook for compiling, running, and plotting
- `mc/` - Dicrectory with Monte Carlo data and C++ code
- `environment.yml` - Environment file for setting up dependencies for Anaconda

In order to run this Notebook, certain python packages need to be installed which is conveniently done
using [Conda](https://www.continuum.io/downloads). From the command line, simply issue

    conda env create -f environment.yml
    source activate cppm
    jupyter-notebook montecarlo.cpp
