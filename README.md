# Dundun2_labels
Code to accompany the paper: 

> Fink, L., Hörster, M., Poeppel, D., Wald-Fuhrmann, E., & Larrouy-Maestri, P. (2026). Reverse-engineering speech and music categorization from a single sound source. *Attention, Perception & Psychophysics*. DOI forthcoming. Preprint: [psyarxiv.com/2635u](psyarxiv.com/2635u)

## Abstract
Classifying whether an auditory signal is music or speech is important for both humans and computational systems. Although
previous literature suggests that music and speech are easily separable categories, common experimental approaches may bias
findings toward this distinction by relying on stimuli from different sound sources and predefined response labels. Here, we
use stimulus material from the dùndún drum – a speech surrogate that can signal either speech-related or musical content. We
first replicate standard speech-music categorization results (N=108). Then, we depart from the typical experimental procedure by asking new participants (N=180) to sort and label the stimulus material, without predefined categories. Hierarchical
clustering of participants’ stimulus groupings reveals multiple organizing dimensions, with the speech–music distinction
reliably present but secondary under label-free conditions. By reverse-engineering the relationship between sorting behavior,
acoustic features, and semantic labels, we characterize how speech–music categorization relates to other salient perceptual
dimensions and how its behavioral prominence depends on task constraints.

**Public significance statement** How people make sense of the sounds around them is a basic question in perception research,
yet many studies assume that sounds naturally fall into categories such as speech or music. This study shows that whether
listeners categorize sounds as speech-like or music-like depends strongly on how they are asked to perform the task, even
when all sounds come from the same physical source. When listeners are not given category labels, other ways of organizing
the sounds become more prominent.

# Citation
If using anything from this repository, please cite the paper above. 

# Code Structure & Dependencies
The core analysis script is in the R Notebook `dundunLabel_analyses.Rmd`. 
This script will produce all figures and statistical analyses found in the actual paper. 
It relies on the packages listed in `dependencies.R` and the custom functions in `functions.r`. As long as you do not change the directory structure of this repository, the *here()* function called in the main script will load all dependencies and the data located in the `data` directory. 

The `acousticAnalyses` folder contains the MATLAB code to extract features from the dundun recordings. It relies on the MDRQA.m helper function (also included). Additionally, the MIR, Audio System, and Signal Processing toolboxes are required.
Stimuli for this project are not publicly available. These scripts are provided for reference only. 

The `supplemental` folder contains analyses added at the request of reviewers and included in the Supplemental Information of the manuscript. One figure, related to the new FAST analysis, is in the main text of the manuscript (Fig. 5). These scripts should be reproducible and rely on the this.path package for setting relative file paths. All fast analyses can be run from the `00_run_all_fast.R` script. Figures S2 and 5 will be produced, alongside relevant data tables. 




