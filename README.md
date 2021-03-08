## A Hard Look At Relativistic Reverberation in MCG--5-23-16 \& SWIFT J2127.4+5654: Testing the Lamp-Post Model

This repository contains the code data associated with the article titled: [A Hard Look At Relativistic Reverberation in MCG--5-23-16 \& SWIFT J2127.4+5654: Testing the Lamp-Post Model](), published in the Astrophysical Journal.

### Abstract
> X-ray reverberation mapping has emerged as a new tool to probe accretion in AGN, providing a potentially powerful probe of accretion at the black hole scale. The lags, along with relativistic spectral signatures are often interpreted in light of the lamp-post model. Focusing specifically on testing the prediction of the relativistic reverberation model, we have targeted several of the brightest Seyfert Galaxies in X-rays with different observing programs. Here, we report the results from two large campaigns with NuSTAR targeting MCG--5-23-16 and SWIFT J2127.4+5654 to test the model predictions in the 3-50 keV band. These are two of three sources that showed indications of a delayed Compton hump in early data. With triple the previously analyzed exposures, we find no evidence for relativistic reverberation in MCG--5-23-16, and the energy-dependent lags are consistent with a log-linear continuum. In SWIFT J2127.4+5654, although a continuum-only model explains the data, the relativistic reverberation model provides a significant improvement to the energy and frequency-dependent lags, but with parameters that are not consistent with the time-averaged spectrum. This adds to mounting evidence showing that the lag data is not consistent with a static lamp-post model. 


### Description
The analysis is organied into several python notebooks, which sometimes call outside functions either from my [toolset package `aztools`](https://github.com/zoghbi-a/aztools) or the `helpers.py` file. The lag calculations are done using the package [`fqlag`](https://github.com/zoghbi-a/fqlag). 

The analysis for each of the two sources is contained in the folders named `mcg5` and `sw2127` respectively. The `helper.py` script is used by both of them.

The analysis of the `xmm` data which is used briefly in some parts of the analysis uses individual scripts located in `{data-folder}/xmm/` (where `data-folder` refers to either `mcg5` or `sw2127`) and `README` file in each case explains what the scripts do

A quick description of the notebooks is as follows:

- `data.md`: Used for data preparation, reduction and the extraction of the spectra and light curves 
- `lc_psd.md`: Used for estimating the power spectra.
- `lag_22l3.md`: Lag calculations using the full 22 energy bins used in section 3.1 in the paper.
- `lag_22l3b.md`: Lag calculations at coarse energy bins (using 11 bins). This first calculates the lag at the published frequencies, then lag-vs-frequency, and lags at the new frequencies


### Data Products:
The data products are available throught the Open Science Framework (links below). Only data that is not available in the heasarc public archive (or can be derived from it in a straightforward manner) is included. 

There are two sets of files, one for each source, named: [`mcg5.tgz` (1GB)](https://osf.io/9m5du/download) and [`sw2127.tgz` (1GB)](https://osf.io/buhd7/download). Each one contains the following folders:
- `nustar/`: contains the folders `{obsid}_p` for each observations. Each folder contains the light curves and spectra from that observation.
- `xmm/`: contains codes and data products for the xmm data, described in details in `xmm/README`
- `spec_xmm_nustar/`: contains the nustar and xmm spectra and spectral fits used to obtain a time-averaged spectral fit.  A discription of the fits is in `spec_xmm_nustar/README`
- `timing/`: contains the results of the psd and lag calculations, including all the modeling with `kynreverb`. The lag spectra are written into pha format where they can be modeled with `xspec`.
    - `timing/*npz` contains the actual lightcurve data used in the calculations.
    - `timing/psd`: contains the lightcurve/psd calculations and modeling.
    - `timing/lag`: contains the lag numbers in `*npz` files, and the lag spectra in the `pha` folder. The plots, including those in the paper are in the `plots` folder.
