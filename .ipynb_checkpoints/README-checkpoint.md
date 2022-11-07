# Mini Thesis Project

## Steps to reproduce results from scratch 
1. (If haven't already) Install miniconda for Python 3.9 or later using instructions from here: https://docs.conda.io/en/latest/miniconda.html
2. Download this repository using `git clone git@github.com:anthonylst6/thesis.git` or clicking Code -> Download ZIP then unzip the folder
3. Open bash shell in home directory of repository then run `conda env create -f env_thesis.yml`
4. (If haven't already) Set up ECMWF CDS API using instructions from here: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-4-DownloadERA5familydatathroughtheCDSAPI
5. Download raw data by entering `conda activate thesis` then running the `data_download.ipynb` notebook in the `scripts` directory
6. (Alternatively) Enter `conda run -n thesis python data_download.py` from the `scripts` directory
7. (For proper rendering of figures) Make sure that a $\LaTeX$ distribution is installed
8. Reproduce results by running the `results.ipynb` notebook in the `scripts` directory for different focus regions by commenting and uncommenting relevant code
9. (Alternatively) Enter `conda run -n thesis python results_wa.py`, `conda run -n thesis python results_ca.py` and `conda run -n thesis python results_sa.py` from the `scripts` directory
10. (If personal computer is limited in RAM) Edit the `results_wa.pbs`, `results_ca.pbs` and `results_sa.pbs` job scripts in the `scripts` directory to be compatible with target HPC facility (the provided scripts were designed for the UNSW Katana computational cluster), then submit these as a batch job

## System requirements
- For Western Australia (WA) raw data and results, up to 25 GB of storage, and up to 60 GB of RAM over 12 hours if using 8 CPUs
- For Central America (CA) raw data and results, up to 25 GB of storage, and up to 60 GB of RAM over 12 hours if using 8 CPUs
- For South America (SA) raw data and results, up to 80 GB of storage, and up to 140 GB of RAM over 12 hours if using 8 CPUs
- The code can also be run with less RAM and fewer CPU cores but this will be slower and will require manual restarting of code everytime RAM limit is reached (the code was designed to pick up from where it left off). But in this case there should at least be 8 GB of RAM for WA and CA results, and at least 24 GB of RAM for SA results.

## High-level description of functions

## Example usage

## Analysing additional ERA5 variables
- (which are not in default)