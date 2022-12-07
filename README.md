# skmocba
Source code for the SK-MOCBA algorithm

Copyright 2020 Sebastian Rojas Gonzalez 

srojas33@gmail.com

sebastian.rojasgonzalez@ugent.be

Source code for the paper: 

%--------------------------------------------------------------------------

Rojas Gonzalez, S., Jalali, H., & Van Nieuwenhuyse, I. (2020). A multiobjective 
stochastic simulation optimization algorithm. European Journal of Operational 
Research, 284(1), 212-226.

%--------------------------------------------------------------------------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the above copyright notice is
retained.

To run any of the algorithms, execute this command in the terminal:

main('-algorithm', @SKMOCBA, '-problem', @ZDT1, '-M', 2, '-D', 5, '-evaluation', 150, '-mode', 1);

where

M: Number of objectives

D: Number of decision variables

evaluation: Number of infill points to sample

mode:

1: Display data

2: Save data


-> This SKMOCBA algorithm code uses partial implementations of PlatEMO to facilitate its usage:

%--------------------------------------------------------------------------

Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
Computational Intelligence Magazine, 2017, 12(4): 73-87".

%--------------------------------------------------------------------------

-> The GLOBAL class encodes the global properties of the experimental setting.

-> The INDIVIDUAL class encodes the individual properties of each design point. 


System requirements:
1) Code has been tested using Matlab 2018b, 2019b, 2020a and 2020b, 
running on 64-bit Linux Debian/Ubuntu, with an Intel i7 VPro CPU with
12 cores in 2.60 GHz and 32 Gb of RAM. 

2) The MOCBA source code is coded in C, based on this reference:

%--------------------------------------------------------------------------

Chen, C. H., & Lee, L. H. (2011). "Stochastic simulation optimization: 
an optimal computing budget allocation (Vol. 1)". World Scientific.

%--------------------------------------------------------------------------

-> C-Matlab bindings are provided to run the code in Matlab via mocba.mexa64

    --> execute mocba(Mat_Obj',Mat_Var'), where
    
        --> Mat_Obj is the matrix of objective values transposed and
        
        --> Mat_Var is the matrix of variance values transposed 

-> These bindings are compiled for 64-bit Linux and haven't been tested on Matlab
 running on Windows or Mac. 


