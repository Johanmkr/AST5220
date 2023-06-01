# An investigation into the power spectrums of matter and radiation

## ABSTRACT
> The cosmic microwave background radiation can be seen across the entire sky, and can provide great insight into the nature of several cosmological phenomena. It is therefore of interest to investigate this theoretically, which is the main aim of this article. We consider linear perturbation to the FLRW cosmology in conformal Newtonian gauge in order to create a pipeline for predicting the power spectra of matter and anisotropies to the CMB. This is done by first calculating the background cosmology, and recombination epoch, including neutrinos. We then ignore neutrinos, polarisation and reionisation when evolving the perturbation equations in time, in order to solve for the power spectra which can be compared with observations. Ignoring these effects yield a significant discrepancy in the result for small scale modes, but provide enough accuracy for the large scale modes. As a result, we obtaine a pipeline through which we can obtain power spectra for different cosmologies and compare them to data, although the analysis is only done for the ΛCDM-cosmology using data from Aghanim et al. (2020).

## Report

The report can be found [here](https://github.com/Johanmkr/AST5220/blob/main/project/tex/cosmology2_report.pdf).

## Structure of report:
### **Background Cosmology**:

According to the cosmological principle, the universe is homogeneous and isotropic on a large scale. Hence, there are no preferred locations nor directions. Furthermore, we may safely assume that the physical laws that govern our local part of the universe is equally valid elsewhere. 

The aim now is to set up the background cosmology, in which the investigation of all interesting phenomena will take place. Setting up the background cosmology is practically equivalent to solving the _Einstein field equation_: $$G_{\mu\nu} = 8\pi GT_{\mu\nu}$$ where $G_{\mu\nu}$ is the Einstein tensor describing the geometry of spacetime, $G$ is the gravitational constant and $T_{\mu\nu}$ is the energy and momentum tensor. This equation relates the geometry and shape of spacetime itself, to its energy content (matter included). There are many solutions to Einstein's field equation, but we want the solution to govern a universe that is spatially isotropic and homogeneous, but may evolve in time. The spacetime metric that satisfies this conditions is the _Friedmann-Lemaître-Robertson-Walker metric_ (FLRW).

We will use this metric in order to describe the background universe, how it may evolve in time, and its history. 

### **Recombination History**: 
The main goal of this section is to investigate the recombination history of the universe. This can be explained as the point in time when photons decouple from the equilibrium of the opaque early universe.  This is known as the _time of last scattering_, and these photons are what we today observe as the CMB. This period of the history of the universe is thus crucial for understanding the CMB. 

We will start by calculating the free _electron fraction_ $X_e$, from which we may find the _optical depth_ $\tau$. This again enables us to compute the _visibility function_, $g$, and the _sound horizon_, $s$. The latter will be of great importance later. 

Recombination happens because the expansion of the Universe cools it down, making the photons less energetic, which in turn make each interaction in the primordial plasma less energetic. At some point, hydrogen atoms are able to form, reducing the number of free electron, hence reducing photon interactions, until they scatter for the last time. We will determine the time of recombination from the free electron fraction, which indirectly tell us how large portion of the free electrons have (re)-combined. Due to the decrease of free electrons, photons interact less with them. At some point, photons scatter for the last time, and this information is encapsulated in the visibility function.

### **Perturbations**:
The aim of this section is to investigate how small fluctuations in the baryon-photon-dark-matter fluid in the early Universe grew into larger structures. This is done by examining the interplay between the fluctuation of these fluids and the subsequent fluctuations of the space-time geometry, originating from tiny quantum fluctuation in the very early Universe. We will model this by perturbing the flat FLRW-metric using the conformal-Newtonian gauge. This will impact how the Boltzmann equations for the different species behave, from which we are able to construct differential equations for key physical observables, and their initial conditions. 


### **Power Spectrums**:
In this section we will construct the angular power spectrum and the matter power spectrum in order to compare our theoretical predictions with actual observables.



## Code explanation

explain

## Run the code

Navigate to the `code`-folder. If in the main directory this is done with:

    cd project/code

### Running c++ code:

The file `Main.cpp` in the `code/src`-folder contain two boolean expression at the top: 

* ``output`` -> if **true**: code should write output files. These are placed in the `code/data` folder as `.csv` files. Default is **true**. 
* ``supernovafit`` -> if **true**: code runs the supernova fit procedure and place results in `code/data` folder.  Default is **false**.

Build c++ codes:

    make cmb

Run:

    ./cmb

Clean:

    make clean

### Plotting
All plotting scripts are located in the `code/plot` folder. The script `plot_utils.py` contain one boolean expression at the top:

* ``TEST``-> if **True**: code shows all the plots rather than saving them as `.pdf`s. Default is **True**. 

Run plots:

    make plots

Clean the plots generated:

    make clean_pdfs




