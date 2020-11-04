# Introduction
This is the custom MATLAB code that accompanies the following publication:

Zuure, M.B., Hinkley, L.B., Tiesinga, P.H.E., Nagarajan, S.S., & Cohen, M.X. (2020). Multiple midfrontal thetas revealed by source separation of simultaneous MEG and EEG. _Journal of Neuroscience, 40_(40), 7702-7713. DOI: https://doi.org/10.1523/JNEUROSCI.0321-20.2020

This code is made available here to comply with the Journal of Neuroscience's code accessibility policy. It is suited only to analyze the data set made available [here](https://doi.org/10.34973/a21y-3t91). I have attempted to make it readable and usable, for verification purposes and for others to learn from. However, it remains academic code and is not fully polished.

This code performs computationally intensive analyses on a large (40 GB) data set. It is recommended that you execute it on a machine with >8 cores and >50 GB of RAM. Some parts, particularly the permutation testing, Granger causality analyses, and noisy eigenvector analysis, may take long (days) to run (but may be sped up, see **How to use**).

## Dependencies
- MATLAB; R2017b was used but earlier versions may work
- EEGlab toolbox, available [here](https://sccn.ucsd.edu/eeglab/index.php). We used v14.1.2.
- Multivariate Granger Causality toolbox v1.0, available [here](http://www.sussex.ac.uk/sackler/mvgc/). We used v1.0.
- bluewhitered plotting package v1.0, available [here](https://nl.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
- distributionPlot package v1.15, available [here](https://nl.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m)
 - the `extractfield` function from the MATLAB Signal Processing toolbox, available from MathWorks [here](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/53545/versions/6/previews/tools/extractfield.m/index.html), to be placed in the `scripts` folder as `extractfield.m`
 - The data set to be analyzed, available [here](https://doi.org/10.34973/a21y-3t91).

## How to use

 1. Clone this repository (or download the files) to a location of your liking. 
 2. Download the data set from [here](https://doi.org/10.34973/a21y-3t91) and put it in the `data` folder, so that the `data` folder contains subfolders `behavior`, `MEGica`, `S01`, etc. 
 3. Change the hardcoded paths in`setpaths.m` to point to the different directories and dependencies listed above.
 4. In MATLAB, navigate to the `scripts` folder and execute the `simon_0x` scripts in order. I recommend doing this in a tmux terminal or similar, as you want to be able to let this run in the background for multiple days; particularly `simon_01_GED_permtest.m`, `simon_04_Granger_causality.m`, and `simon_06_eigennoise.m` will be slow. Note: for faster but less accurate significance thresholding, change the `num_perm` parameter in `simon_01_GED_permtest.m`from 1000 to something lower. This will affect the number of significant components found per subject. To speed up `simon_06`, decrease the `iterations` parameter from 50 to something lower (decreasing the number of correlations computed per noise level), or consider running multiple instances of MATLAB in parallel, each executing the code for a single subject.
 5. After everything has been run, the `figures` folder contains the .pdfs that the manuscript figures are created from. The intermediate generalized eigendecomposition results and further analysis results are in the `results` folder, in the `S0x_GED.mat` and `S0x_ana.mat`files, respectively. These can be loaded into MATLAB and inspected. For context on what these (many) variables are, please refer to the code that created them; I've done my best to keep this readable.

## Summary of scripts

* `simon_00_behavioral` applies a two-way ANOVA to trial reaction times and accuracy, to identify congruency sequence effects.
* `simon_01_GED_permtest` extracts the M/EEG data features of interest and enters them into a generalized eigendecomposition (GED), returning eigenvalues and eigenvectors (components). If no previous permutation test results exist, the GED is repeated an additional 1000 times per subject with shuffled data to determine eigenvalue significance---this is a slow process! For faster but less accurate significance thresholding, change the `num_perm` parameter in this file from 1000 to something lower, maybe 100 to 300 for passable accuracy.
* `simon_02_explore_midfrontal` takes the significant components identified by permutation testing, selects for midfrontalness, and extracts component features such as theta power over time and between-component synchrony.
* `simon_03_sensors` performs sensor-level analyses (as opposed to component-level analyses), such as time-frequency decomposition.
* `simon_04_Granger_causality` computes Granger causality over time between each subject's component time series, estimating information transfer over time between components.
* `simon_05_PCA` performs PCAs on the theta power time series, GC time series, and MEG topographies for all components, extracting common temporal and spatial features.
* `simon_06_eigennoise` creates noisy eigenvector copies and computes the pairwise correlations, to test the ground truth correlations against noise-level-specific distributions.
* `simon_07_publication_plots` generates and saves the figure panels used in the manuscript. Resulting PDFs were optimized for editing in Adobe Illustrator and may appear small or cropped, or be split into multiple parts.
* `simon_08_statistics` performs additional statistical tests and ad hoc inspections, some of which appear in the manuscript. This script mainly produces terminal output and renders some figures, and is best executed cell by cell in MATLAB.
* `r2r_ICA`, `r2r_split_half_validity`, and `r2r_multimodality` perform robustness analyses requested by the reviewers. These scripts mainly produce terminal output and render some figures, and are best executed cell by cell in MATLAB.

## License
This code is released under the MIT License.

## Contact
For further questions regarding this code, please refer to the manuscript or contact me (Marrit Zuure) at my github-associated e-mail address.
