# Unlabeled Compressed Sensing (UCS)
##### Mohamed Akrout, Faouzi Bellili, Amine Mezghani
This repository contains the Matlab code of the unlabeled compressed sensing (UCS) algorithm proposed in the paper: [Unlabeled Compressed Sensing from Multiple Measurement Vectors](http://arxiv.org/abs/2406.08290)
## Abstract
This paper tackles a general form of the unlabeled compressed sensing (UCS) problem with multiple measurement vectors (MMV). The goal is to recover an unknown structured signal matrix, $\mathbf{X}$, from its noisy linear observation matrix, $\mathbf{Y}$, whose rows are further randomly shuffled by an unknown permutation matrix $\mathbf{U}$. A new Bayes-optimal UCS recovery algorithm is developed from the bilinear approximate message passing (Bi-VAMP) framework using non-separable and coupled priors on the rows and columns of the permutation matrix $\mathbf{U}$. In particular, standard unlabeled sensing is a special case of the proposed framework, and UCS further generalizes it by neither assuming a partially shuffled signal matrix $\mathbf{X}$ nor a small-sized permutation matrix $\mathbf{U}$. For the sake of theoretical performance prediction, we also conduct a state evolution (SE) analysis of the proposed algorithm and show its consistency with the asymptotic empirical mean-squared error (MSE). Numerical results demonstrate the effectiveness of the proposed UCS algorithm and its advantage over state-of-the-art baseline approaches in various applications. We also numerically examine the phase transition diagrams of UCS, thereby characterizing the detectability region as a function of the signal-to-noise ratio (SNR).

## Repository Structure
This repository contains two folders:
  - **UCS**: it contains the code of the UCS algorithm solving the unlabeled sensing problem whose observation model is $\boldsymbol{Y}~ =~ \boldsymbol{U}\boldsymbol{A}\boldsymbol{X} + \boldsymbol{W}$. This folder has two files:
    * *UCS.m*: the file containing the algorithmic steps of UCS.
    * *UCS_opt.m*: the file defining the parameters of the UCS simulation.
  - **priors**: it contains the code of the prior modules employed by the UCS algorithm. This folder has three files:
    * *column_prior_U.m*: the column prior on $\mathbf{U}$.
    * *row_prior_U.m*: the row-wise assignment prior on $\mathbf{U}$.
    * *prior_gauss_X.m*: the Gaussian prior on $\mathbf{X}$.

## Running Experiments
- To run UCS, execute the file `main.m`. The template of the output in the Matlab console is:
```
===== Problem dimension =====
 - N = 100
 - M = 50
 - R = 10

===== Running UCS =====
[t=0] nrmse =  
[t=100] nrmse =
[t=200] nrmse = 
.
.
.
Final nrmse = 
Running time = 
```
## Important facts
- We use a damping factor defined in the [UCS_opt](https://github.com/makrout/Unlabeled-Compressed-Sensing/blob/main/UCS/UCS_opt.m) object and a variance annealing parameter defined [here](https://github.com/makrout/Unlabeled-Compressed-Sensing/blob/main/main.m#L19). Both of them should be tuned depending on the problem dimension (i.e., n, m, and r).
- We are impressed how UCS can find a permutation matrix by approximating the permutation prior using a cascade of two assignment priors. The simulation in the file `main.m` finds a permutation matrix of size $N \times N$ for $N=50$.
- The results of the UCS algorithm in the paper are reported using Monte Carlo simulations. This means that when the UCS algorithms fails to converge sometimes in one specific simulation, it does not mean that the UCS algorithm does not work. You can run it multiple times, change the initialization values of variables, the damping factor, and the variance annealing parameter.
- *Why is Bi-VAMP the adequate framework for UCS*: Unlike all existing competitors for bilinear recovery (e.g., BAd-VAMP, Bi-GAMP, and P-BiG-AMP), Bi-VAMP is the only bilinear algorithm that is based entirely on a factor graph consisting of vector-valued variable nodes in *both* $\mathbf{U}$ and $\mathbf{V}$. This unique feature of BiG-VAMP makes it possible to use two *inseparable* priors on the columns and rows of $\mathbf{U}$ so as to enforce the overall permutation prior. In other words, the factor graphs of the existing bilinear algorithms consist of scalar-valued variable nodes and hence do not capture the correlation inside the rows/columns of $\mathbf{U}$.
  
## Citing the paper (bib)

If you make use of our code, please make sure to cite our paper:
```
@misc{akrout2024unlabeled,
      title={Unlabeled Compressed Sensing from Multiple Measurement Vectors}, 
      author={Mohamed Akrout and Amine Mezghani and Faouzi Bellili},
      year={2024},
      eprint={2406.08290},
      archivePrefix={arXiv},
      primaryClass={cs.IT}
}
```
