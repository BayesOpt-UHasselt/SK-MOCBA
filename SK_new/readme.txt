
Disclaimer: 

This software has been edited and extended by Sebastian Rojas Gonzalez 
<sebastian.rojasgonzalez@ugent.be><srojas33@gmail.com>, based on on the source code provided by 

Ankenman, B., Nelson, B. L., & Staum, J. (2010). 'Stochastic kriging for simulation metamodeling.'
Operations research, 58(2), 371-382.

The original software is copyrighted by Barry L.Nelson, Jeremy Staum, Evern Baysal and Wei Xie 2009 (see
original reference for further copyright details). 

========================= Instructions ======================

The purpose of SK software is to provide an accurate estimate of the response surface 
of interest. Given the simulation outputs at design points, this software can be used 
to compute the Maximum Likelihood Estimation (MLE) of model parameters through SKfit_new.m, 
and to provide an accurate estimate of response surface at prediction points through 
SKpredict_new.m. 

In the following:
k represents the number of design points
K represents the number of prediction points
d represents the dimensionality of  design variables
b represents the number of basis functions for trend estimation.


1. SKfit_new.m
This function can be used to fit a stochastic kriging model to simulation outputs at design 
points, and extract MLE of model parameters. 

	model = SKfit_new(X, Y, B, Vhat, gammaP,algor_sel)

Inputs: include matrix of design points X (k x d), simulation outputs at design points 
(mean Y (k x 1) and intrinsic variance Vhat (k x k)), matrix of basis functions B (k x b) for trend 
estimation, and type of correlation function (gammaP = 1: exponential, 2: gauss, 3: cubic, 4: Matern 2.5). 
The optimization can be done by the following methods (algor_sel = 1: trust-region-reflective algorithm,
2: active-set algorithm, 3: interior-point algorithm).

Matrix B is the design matrix used to approximate the trend of response surface. 
For example, for cubic trend, the i-th row of B is [1 xi xi^2 xi^3], where xi is the i-th design point. 
Without any preliminary information about trend, the default B is a k-length column unit vector.
		 
Outputs: model gives the MLEs of model parameters, including tausquared, 
theta, beta as shown Eq.(13) in the reference, and some intermediate variables which will be used 
for SK prediction. 
	

2. SKpredict_new.m
This function can be used to make response surface estimations at prediction points based on 
MLE parameter estimates obtained from SKfit_new.m.

	f = SKpredict_new(model,Xpred,Bpred)

Inputs: include the output from SKfit_new.m model, matrix of prediction points Xpred (K x d), 
matrix of basis function for trend estimations at prediction points Bpred (K x b). 
For example, for the cubic trend, the i-th row of Bpred is [1 xi xi^2 xi^3], 
where xi is the i-th prediction point. 

Output: f represents the fitted values at prediction points (K x 1), Eq.(12) in reference


Reference: 

Jeremy Staum, "Better Simulation Metamodeling: the Why, What, and How of Stochastic Kriging", 
2009 Winter Simulation Conference.


====================== Troubleshooting =========================

SK software is designed to track two kinds of errors, the mismatch of function inputs 
and singularity of covariance matrix. When the smooth correlation function is chosen, 
and the intrinsic variance is very small, singularity problem may show up as 
the number of design points increases. In the SK software, matlab function fmincon.m 
is used to iteratively search for MLE of model parameters that includes taking inverse of 
covariance matrix of outputs. When the covariance matrix is near singular, the estimated 
covariance matrix may not be positive definite. Following error message will show up:

	Error using ==> logPL
	covariance matrix is nearly singular

Note because of the singularity issue, at this stage this software is not recommended
to use to fit kriging model with smooth correlation functions. In addition, when the 
number of design points becomes very large and the dimensionality of the design variables 
increasing, SK software may need much longer computation time.


======================== Example ==============================
M/M/1 Queue 

MM1sim.m - simulate M/M/1 queueing problem
MM1_example.m - Based on M/M/1 simulation outputs, fit SK model and make predictions

MM1_example.m is used to fit the response surface of expected waiting time in M/M/1 
queue over service rate. The arrival rate is fixed as 1. Design variable is system utilization.

This example uses stochastic kriging with gauss correlation function and constant trend estimation. 
The first part (lines 12~32) is to generate the evenly distributed design and prediction points, 
and obtain the simulation outputs at design points. Then, given the simulation outputs 
(Y,Vhat), "SKfit_new" is used to obtain the parameter estimates. After that, 
based on model parameter estimates, fitted values at prediction points can be computed by 
"SKpredict_new".

