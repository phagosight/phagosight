# Phagosight
Tracking Algorithms for Phagocytes and other moving objects.
## Installation
After cloning the repository, it is important to add it (with subfolders) to the
Matlab path. Full manual and instructions can be found
[here](http://www.phagosight.org/).

## Usage and Quick Start
[Phagosight](http://www.phagosight.org/) can be used with the help of a GUI that
allows the user to select manually the folder where the data is in. This is done
by the command
```Matlab
[handles] = neutrophilAnalysis();
```
The program can be fully automatic, the user knows some or all the parameters.
An example of this:
```Matlab
[handles] = neutrophilAnalysis('/path/to/data/', 0, [], [0.2 0.4]);
```
where parameters not known, or that want to be taken by default, are simply
set to empty arrays `[]`. For a detailed review of the input parameters, the
`help` command outputs all the information the user might need.
```Matlab
help neutrophilAnalysis
```
### User Manual

For a detailed and comprehensive user manual, check the Wiki section.

### Output and results
Phagosight returns variable `handles`, which contains information about the 
tracks and locations where results from intermediate steps are stored. For 
example, if the path to data is a folder with name `Data/`, then the new 
folders created:
* `Data_mat_Or` Stores the original data in Matlab's own `.mat` format.
* `Data_mat_Re` Stores the data after preprocessing, whether it is reduction 
 by subsampling or filtering (if size is not reduced)
* `Data_mat_La` Stores the labelled data after segmentation.
## Citation
Phagosight is released under the GNU General Public License v3. 
Please cite Phagosight in your publications if it helps your research:
```BibTex
@article{Henry2013,
 author = {Henry, Katherine and Pase, Luke and Ramos-Lopez, Carlos Fernando and
 Lieschke, Graham J. and Renshaw, Stephen A. and Reyes-Aldasoro, Constantino Carlos},
 issn = {19326203},
 journal = {PLOS ONE},
 month = {January},
 number = {8},
 pages = {e72636},
 pmid = {24023630},
 title = {{PhagoSight: An Open-Source MATLAB\textsuperscript{\textregistered}\
 Package for the Analysis of
 Fluorescent Neutrophil and Macrophage Migration in a Zebrafish Model}},
 volume = {8},
 year = {2013}
}
```
 
