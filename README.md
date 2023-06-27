# Dundun2_labels
Code to accompanying the paper: 

> Fink, L., Hörster, M., Poeppel, D., Wald-Fuhrmann, M., Larrouy-Maestri, P. (2023, submitted). Features underlying speech versus music as categories of auditory experience. 

If using anything from this repository, please cite the paper. 

# Contents

The core analysis script is `dundunLabel_analyses.Rmd`. This script will produce all figures and statistical analyses found in the paper. It relies on the packages listed in `dependencies.R` and the custom functions in `functions.r`. As long as you do not change the directory structure of this repository, the here() function called in the main script will load all dependencies and the data located in the `data` directory. 

The `acousticAnalyses` folder contains the MATLAB code to extract features from the dundun recordings. It relies on the the MDRQA.m helper function (also included). Additionally, the MIR, Audio System, and Signal Processing toolboxes are required.

[All stimuli are available for download here](https://edmond.mpdl.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.JATDRF)
