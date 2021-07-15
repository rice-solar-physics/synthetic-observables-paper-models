# Understanding Heating in Active Region Cores through Machine Learning I. Numerical Modeling and Predicted Observables

- [arXiv](https://arxiv.org/abs/1906.03350)
- [ApJ](https://iopscience.iop.org/article/10.3847/1538-4357/ab290c)
- [ADS](https://ui.adsabs.harvard.edu/abs/2019ApJ...880...56B)

## Authors

* Will Barnes, *Department of Physics & Astronomy, Rice University*
* Stephen Bradshaw, *Department of Physics & Astronomy, Rice University*
* Nicholeen Viall, *NASA Goddard Space Flight Center*

## Abstract

To adequately constrain the frequency of energy deposition in active region cores in the solar corona, systematic comparisons between detailed models and observational data are needed. In this paper, we describe a pipeline for forward modeling active region emission using magnetic field extrapolations and field-aligned hydrodynamic models. We use this pipeline to predict time-dependent emission from active region NOAA 1158 as observed by SDO/AIA for low-, intermediate-, and high-frequency nanoflares. In each pixel of our predicted multi-wavelength, time-dependent images, we compute two commonly-used diagnostics: the emission measure slope and the time lag. We find that signatures of the heating frequency persist in both of these diagnostics. In particular, our results show that the distribution of emission measure slopes narrows and the mean decreases with decreasing heating frequency and that the range of emission measure slopes is consistent with past observational and modeling work. Furthermore, we find that the time lag becomes increasingly spatially coherent with decreasing heating frequency while the distribution of time lags across the whole active region becomes more broad with increasing heating frequency. In a follow up paper, we train a random forest classifier on these predicted diagnostics and use this model to classify real AIA observations of NOAA 1158 in terms of the underlying heating frequency.

## Building the Paper

This repository contains all of the code and data needed to reproduce all of the figures in this paper, including training the random forest model that is used for making the heating frequency predictions.

First, build the conda environment and install all of the needed packages. This assumes you're using the [Anaconda distribution of Python](https://www.anaconda.com/products/individual), but this can be done using any distribution or environment manager.

```shell
$ conda create --name barnes-bradshaw-viall-20 --file conda-env.txt
$ conda activate barnes-bradshaw-viall-20
```

Next, build the paper using LaTeX and PythonTeX. You will need to install the PythonTeX package. This will run all of the needed Python code (including inline in the LaTeX source) to build all of the figures and train all of the models.

```shell
$ cd paper
$ pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
$ pythontex --interpreter python:python paper.tex
$ bibtex paper
$ pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
$ pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
```

This may take several minutes. You can find the built PDF in `paper/paper.pdf`.