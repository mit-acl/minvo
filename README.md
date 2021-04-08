# MINVO Basis: Finding Simplexes with Minimum Volume Enclosing Polynomial Curves #

Example for n=3:

[![](./imgs/minvo3d.png)](https://www.youtube.com/watch?v=f_JOYud9LUU) 

## Citation

When using MINVO, please cite this paper ([pdf](https://arxiv.org/abs/2010.10726), [video](https://youtu.be/x5ORkDCe4O0)):

```bibtex
@article{tordesillas2020minvo,
  title={{MINVO} basis: Finding simplexes with minimum volume enclosing polynomial curves},
  author={Tordesillas, Jesus and How, Jonathan P},
  journal={arXiv preprint arXiv:2010.10726},
  year={2020}
}
```

## Instructions to use the MINVO basis

* **For a particular curve (Problem 1 of the paper)**: See an example (with `n=3`) in [`curve_given3D.m`](https://github.com/mit-acl/minvo/blob/master/src/curve_given3D.m)

* **For a particular simplex  (Problem 2 of the paper)**: See an example (with `n=3`) in [`simplex_given3D.m`](https://github.com/mit-acl/minvo/blob/master/src/simplex_given3D.m)

To obtain all the figures of the paper, you can  simply run [`plot_solution.m`](https://github.com/mit-acl/minvo/blob/master/src/plot_solution.m) 

You don't need to install any external solvers for the steps above.

## Instructions to derive the MINVO basis

* The file [`general_formulation_sos.m`](https://github.com/mit-acl/minvo/blob/master/src/general_formulation_sos.m) solves Problem 3 in the paper using SOS.
* The file [`general_formulation_lukacs_theorem.m`](https://github.com/mit-acl/minvo/blob/master/src/general_formulation_lukacs_theorem.m) solves Problem 3 in the paper using the Markov–Lukács Theorem. 
* The file [`formula_formulation.m`](https://github.com/mit-acl/minvo/blob/master/src/formula_formulation.m) solves Problem 4 in the paper. 

Follow the instructions in each file to prove local/global optimality.

Depending on the settings you choose in each file, you may need to install [YALMIP](https://yalmip.github.io/) (tested with [this release](https://github.com/yalmip/YALMIP/releases/tag/R20200116_hotfix)), [SNOPT](https://ccom.ucsd.edu/~optimizers/) and/or [MOSEK](https://www.mosek.com/). 


---------

> **Approval for release**: This code was approved for release by The Boeing Company in December 2020. 
