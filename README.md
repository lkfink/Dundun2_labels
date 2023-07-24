# Dundun2_labels
Code to accompany the paper: 

> Fink, L., Hörster, M., Poeppel, D., Wald-Fuhrmann, M., & Larrouy-Maestri, P. (2023, July 24). Features underlying speech versus music as categories of auditory experience. Retrieved from [psyarxiv.com/2635u](psyarxiv.com/2635u)

## Abstract
Listeners show remarkable abilities ​to distinguish music and speech when asked but the essence of such categories is arguable. Here, using recordings of dùndún drumming (a West-African drum also used as a speech surrogate), we first replicate standard speech-music categorization results (N=108, sample size based on a prior study), then depart from the typical experimental procedure by asking participants (N=180) to freely categorize and label these recordings. Hierarchical clustering of participants’ stimulus groupings shows that the speech/music distinction emerges, but is not primary. Analysis of participants' labels in the free-response task converges with acoustic predictors of the categories, supporting the effect of priming in music/speech discrimination, and thereby providing a new perspective on the categorisation of such common auditory signals. 

# Citation
If using anything from this repository, please cite the paper above. 

# Code Structure & Dependencies
The core analysis script is in the R Notebook `dundunLabel_analyses.Rmd`. 
This script will produce all figures and statistical analyses found in the paper. 
It relies on the packages listed in `dependencies.R` and the custom functions in `functions.r`. As long as you do not change the directory structure of this repository, the *here()* function called in the main script will load all dependencies and the data located in the `data` directory. 

The `acousticAnalyses` folder contains the MATLAB code to extract features from the dundun recordings. It relies on the MDRQA.m helper function (also included). Additionally, the MIR, Audio System, and Signal Processing toolboxes are required.

All stimuli [are available for download here.](https://edmond.mpdl.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.JATDRF)




