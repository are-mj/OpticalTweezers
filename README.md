# OpticalTweezers
Matlab functions for post-processing of protein stretching experiments

This repository contains Matlab code to analyse experiment files from optical tweezer experiments stretching and unfolding single protein molecules, followed by relaxing the force until refolding.
This is an improved version of the one in the repository "protein stretching". 

Central files:
analyse_experiment.m: Reads an experiment file and automatically identifies an calulates properties of rips and zips
analyse_many.m:  Run analyse_experiment for a list of files
RipAlanysis.mlapp:    App for inspecting and modifying results from analyse_experiments

The results are given as Matab tables Tp (pulling traces) and Tr (relaxing traces)
