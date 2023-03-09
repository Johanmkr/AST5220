# AST5220 - Cosmology II - Numerical project

## The work done in order to calculate the CMB power spectrum. 

## Report

The report can be found [here](https://github.com/Johanmkr/AST5220/blob/main/project/tex/cosmology2_report.pdf).

## Milestone 1 - Background cosmology

In milestone 1 we set up the background cosmology and solve for and isotropic and homogenous universe. 

## Milestone 2
some

## Milestone 3
some

## Milestone 4
some



## Code explanation

explain

## Run the code

Navigate to the `code`-folder. 

Comment/uncomment line 52-54 in `main.cpp` in order to disable/enable the supernova-fitting.

Build c++ codes:

    make cmb

Run:

    ./cmb

Clean:

    make clean

Run plots for testing:

    make plot_test

Run plots for analysis:

    make plot_analM1

Clean the plots generated:

    make clean_pdfs



## TODO:
 *  Optimise code
 *  Derive Friedmann equations in an appendix.


 ### Milestone1
 *  Plots
    * Show analytical values in different regimes (test plots)
    * Conformal Hubble, show vertical lines and acc. onset. 
    * More sensible units of conformal Hubble.
    * COSMIC TIME -> log yaxis.
    * ERAS -> use vertical lines
    * Times -> update table to be more readable.
    * Give value of 1$\sigma$ constraint.
    * Use "fiducial cosmology" instead of prediction.

* Code
    * Use OmegaM and OmegaR_tot

* Report
    * Write introduction
    * Write Abstract
    * Update nomenclature

* Theory
    * Use "density" instead of "density of curvatere" (only mathematical similaritites, not a real density)
    * Redshift - what is it really?
    * Conformal time $\eta$ - what is it really, why do we care about it?
    * Radial distance is really the radial coordiante in the FRLW line element. 
    * Typo in eq. 19
    * Derive acceleration onset. 

* Sanity checks
    * Nice, use \ll as $\ll$ instead of << 

* Methods
    * Some sections should be merged and renamed. ODEs and spline are good. 
    * Naming: model evaluation -> something, cut the resemble reality stuff.
    * Fix typos.
    * $\chi^2$ not really errors, but tells us something of the likelihood. 

* Results
    * Talk more about $\sigma$ intervals.
    * Quote constraints on different parameters.
    * Discuss physical interpretation.

* Appendix
    * Check feedback on devilry and fix when time is available (lol)


### Milestone 2
* Start


 ## Questions
 * FIXED -> no use more sensible. Right unit on the $\mathcal{H}$ plot?
 * FIXED -> Include lines insteadRegime shading on all plots?
 * FIXED -> appears so. Colours on plots okay?
 * FIXED -> Yay. Nomenclature in that style yay/nay?



