# `BFMAN`

## Introduction

`BFMAN` is a repository that implements sparse Bayesian factor model with mass-nonlocal factor scores. It specifies shrinkage prior on factor loadings and sparse prior on factor scores. These settings additionally induce sparsity, allowing more flexible structure on latent factors. Compared to other factor models, BFMAN demonstrates better performance in estimating loadings and scores, showing more accuracy in factor selection.

## File Structure

`/R code/src` contains required functions to implement BFMAN:  
	&emsp;`../func_lib.R` includes `pmom()` to obtain the density of pMOM distribution, `pos_eta()` to obtain the posterior of scores $\eta$, and `MH()` as a sampling algorithm.  
	&emsp;`../mdl.R` includes codes of BFMAN and MGPS (multiplicative Gamma process shrinkage prior) for parallel jobs. It is called by our simulation study.  

`/R code/simul` contains files for simulation.  
	&emsp;under `../scripts` folder:  
		&emsp;&emsp;`../data_gen.R` is used for data generation. For specific settings, please refer to our paper.  
		&emsp;&emsp;`../main.R` is used for calling different models with replications.  
		&emsp;&emsp;`../simul_debug.R` is a simplified script for running BFMAN.  
		&emsp;&emsp;`../eval_perform.R` is used to reorganize the outputs and compare the performance of various models.  
	&emsp;`../data` folder contains simulated data with different scenarios.  
	&emsp;`../results` folder stores output of BFMAN and MGPS.  

`/R code/nutrients` contains files for real data analysis.
	&emsp;under `../scripts` folder:
		&emsp;&emsp;`../nut_analy.R`  is for applying BFMAN on HCHS/SOL data set.  
		&emsp;&emsp;`../nut_vis.R` is for visualizing output of BFMAN, including heatmap of loading matrix $\Lambda$, and OR estimation via Bayesian generalized linear regression.  
	&emsp;`../data` folder contains the data files that includes nutrient consumptions, health records and demographic information.  
	&emsp;`../results` folder stores output of `nut_analy.R`.  

## Quick start

`simul_debug.R` is a simple script to implement BFMAN. 

After running the full script, `Lambda_hat` is the estimated $\Lambda$ obtained with median estimator, so does `eta_hat`. We use RV coefficient to evaluate the performance of estimation.

## Model Specification

For generic form of factor model:
$$
y_i = \Lambda \eta_i + \epsilon_i, (i = 1,\cdots, n)
$$
$y_i \in \mathbb{R}^p$ observed variables
$\Lambda \in \mathbb{R}^{p \times k}$: factor loadings
$\eta_i \in \mathbb{R}^k$: latent factor scores
$\epsilon_i \sim \mathcal{N}_p(0, \Sigma)$ is the error term with covariance $\Sigma = \mbox{diag}(\sigma_1^2, \cdots, \sigma_p^2)$

we specify MGPS on factor loadings $\Lambda$:
$$
\lambda_{jh} \mid \phi_{jh}, \tau_h \sim \mathcal{N}(0, \phi_{jh}^{-1} \tau_h^{-1}) \\
\phi_{jh} \sim \mathrm{Ga}(\nu/2, \nu/2)\\
\tau_h = \prod_{l=1}^{h} \delta_l, \quad
\delta_1 \sim \mathrm{Ga}(a_1, 1), \quad 
\delta_l \sim \mathrm{Ga}(a_2, 1), \quad l \ge 2
$$
Mass-nonlocal prior on scores $\eta_i$:
$$
\eta_{ih} \sim (1 - Z_{ih}) \delta_0(\eta_{ih}) + Z_{ih} \ \mbox{pMOM}(\eta_{ih})
$$
$$
Z_{ih} \sim \mbox{Bern}(\theta_h), \theta_h \sim \mbox{Beta}(\alpha, \beta)
$$

A Gibbs sampler embedded with Metropolis-Hastings algorithm is derived for parameter estimation.
