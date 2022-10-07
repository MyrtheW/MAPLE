# MAPLE sequencing errors
Useful files: 
- `MAPLEv0.1.9_error_site_specific.py` This file is capable of inferring topologies without taking into account error rates, with error rates or with site specific error rates.
  This file incorporates all updates from v0.1.9 that were marked by `TODO Myrthe`. Important arguments are: '--input', '--output', '--errorRate' and '--siteSpecific'. Please refer to there help descriptions for more information.
- `MAPLEv0.1.9_error_site_specific_merged.py` I created this file by comparing the file `MAPLEv0.1.9` of the 3rd of October and merged all other changes into my error version. 
  However, when trying to call the script with an error rate, a `None` vector is created and causing issues that still need to be debugged. I stopped with this for now, as you are still working on version 0.1.9.
- `MAPLE_simulate_errors.py` This script can be used to simulate (position specific) errors. Important arguments are: '--input', '--output', '--errorRate' and '--siteSpecific'. Please refer to there help descriptions for more information. 
- `benchmarking_bash_scripts.py` This file creates bash scripts to create fasta files with (possibly position specific) simulated errors, create a MAPLE file from those and call my version of MAPLE. File names, number of repeats and sample size should be easy to adjust. When a benchmarking file is provided, my MAPLE version stores the relevant results in here. 
- `benchmarking_result_analysis.py` This is (not a very neat) python script used for making the plots from the benchmarking results (the tsv file).
