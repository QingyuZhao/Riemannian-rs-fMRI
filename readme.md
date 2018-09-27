# A Riemannian Framework for Longitudinal Analysis of Resting-State Functional Connectivity

This page gives a toy example of performing longitudinal analysis on rs-fMRI connectivity based on 
the Riemannian geometry of covariance matrices. In this scenario, there are two groups of subjects. Each subject has multiple covariance matrices assumed to be estimated from that subject's longitudinal rs-fMRI data. The task is to perform group analysis on the longitudinal trajectory of covariance matrices, for each matrix entry. 

See our MICCAI paper

A Riemannian Framework for Longitudinal Analysis of Resting-State Functional Connectivity, Q Zhao, D Kwon, KM Pohl, International Conference on Medical Image Computing and Computer-Assisted Intervention 2018

for model description and experimental explaination. 

# Run the code
I ran these codes using Matlab 2016b

RUN demoExample.m to perform

1. Computing Subject-Specific Longitudinal Trajectories

2. Group Analysis on the trajectories on each matrix element


