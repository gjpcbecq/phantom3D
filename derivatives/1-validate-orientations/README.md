# Resources and code for validating the 3D phantom orientation

<a href="./figs/patterns.png"><img src="./figs/patterns.png" width=300></a>

## Introduction

This repository provides code and resources associated with the results presented in the article _A 3D-printed phantom to validate subject orientation in 3D imaging and recordings_. 
The resources allow you to:
1. __Generate__ the 3D phantom template for CT and MRI recordings.
2. __Get the best scores__ for flips and permutations on analyzed volumes.
3. __Correct__ the affine transform in incorrect NIfTI recordings.


Part of the analyzed dataset is available on OpenNeuro at https://doi.org/10.18112/openneuro.ds006123.v1.0.2

## Requirements 

The code and results were generated using the following software versions: 
- python 3.12
- numpy 2.4.0
- matplotlib 3.10.8
- nibabel 5.3.3
- ImageIO (imageio.v3) 2.37.2
- pandas 2.3.3
- jupyterlab 4.5.1

## Content

- README.txt` – this file
- `code/` –  Code to generate resources in `figs/` and `templates/`.
- `templates/` – NIfTI files for templates.
- `figs/` – Content to generate figures for the article.

## List of files in `code/` 

### Python modules 
- `tools.py` – Custom functions developped to perform analyses. 
- `ext.py` – Imports external modules so that functions in `tools.py` are prefixed with `ext`. 

### Jupyterlab notebooks 

__Note__: The database pathname in the Python notebooks must be updated to reproduce the results. 
Processing times may vary depending on your system configuration. 

- `1-compute_patterns-ct.ipynb` – Generate CT pattern       
- `1-compute_patterns-mri.ipynb` – Generate MRI patterns
- `2-check_source_sub-01_ses-1_gin.ipynb` – Best fit and corrections for Table 2 (1a). 
- `2-check_source_sub-01_ses-2.ipynb` – Best fit + corrections for Table 2 (2).
- `2-check_source_sub-01_ses-4_CORO.ipynb` – Best fit + corrections for Table 2 (4).
- `2-check_source_sub-01_ses-5_CORO.ipynb` – Best fit for Table 2 (5). 
- `3-check_sub-01_ses-1_brkraw.ipynb` – Best fit for Table 2 (1b).
- `3-check_sub-01_ses-1.ipynb` – Validates `sub-01_ses-1`.
- `3-check_sub-01_ses-2.ipynb` – Validates `sub-01_ses-2`.
- `3-check_sub-01_ses-4.ipynb` – Validates `sub-01_ses-4`.
- `3-check_sub-01_ses-5.ipynb` – Validates `sub-01_ses-5`.
- `4-figures.ipynb` – Generate figures for documentation and the article.
- `4-compile_tex.ipynb` – Compile the TeX files.
- `4-example.ipynb` – Examples for `tools.py`.

### LaTeX Files for Figure Generation (Using TikZ)

- `4-fig4.tex` – Generates `fig4.pdf`
- `4-fig5.tex` – Generates `fig5.pdf`
- `4-fig4.pdf`
- `4-fig5.pdf`
- `4-figures.tex` – TikZ commands to generate Fig.4 and Fig.5.  

### Others
- `LICENSE.txt` – The detailed license agreement in use for all code in `code` folder.  
- `/img` – Contains masks images, and the GIMP project used to generate masks. 


__Author__: GJC Becq  
__Date__: 2026-04-03  
__Copyright__: CeCILL (compatible GNU see `LICENSE.txt`)  



