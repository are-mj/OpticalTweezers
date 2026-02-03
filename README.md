# OpticalTweezers
Matlab functions for post-processing of protein stretching experiments

This repository contains Matlab code to analyse experiment files from optical tweezer experiments where a single protein is stretched until unfolding (rip) and then allowed to relax and refold (zip)  

Central files:

analyse_experiment.m: Reads an experiment file and automatically identifies an calulates properties of rips and zips

analyse_many.m:  Run analyse_experiment for a list of files

RipAlanysis.mlapp:    App for inspecting and modifying results from analyse_experiments

The results are given as Matab tables Trip (pulling traces) and Tzip (relaxing traces)
