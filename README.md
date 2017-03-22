# phagosight
Tracking Algorithms for Phagocytes and other moving objects.

## Instalation
After cloning the repository, it is important to add it (with subfolders) to the
Matlab path. Full manual and instructions can be found
[here](http://www.phagosight.org/).

## Usage and quick start
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
