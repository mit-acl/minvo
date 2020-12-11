# MINVO Basis: Finding Simplexes with Minimum Volume Enclosing Polynomial Curves #

Example for n=3:

![](./imgs/minvo3d.png) 

## Citation

When using MINVO, please cite [this paper](https://www.google.com/):

```bibtex
@article{tordesillas2020minvo,
  title={{MINVO} Basis: Finding Simplexes with Minimum Volume Enclosing Polynomial Curves},
  author={Tordesillas, Jesus and How, Jonathan P},
  journal={arXiv preprint},
  year={2020}
}
```

## Instructions to use the MINVO basis

* **For a particular curve (Problem 1 of the paper)**: See an example (with `n=3`) in [`curve_given3D.m`](https://github.com/mit-acl/minvo/blob/master/src/curve_given3D.m)

* **For a particular simplex  (Problem 2 of the paper)**: See an example (with `n=3`) in [`simplex_given3D.m`](https://github.com/mit-acl/minvo/blob/master/src/simplex_given3D.m)


## Instructions to derive the MINVO basis

#### To generate all the figures of the paper
You can simply run `plot_solution.m` (you don't need to install any external solver to do this).

#### To run the optimization
* The file `general_formulation.m` solves Problem 3 in the paper. Follow the instructions in that file to prove local/global optimality.
* The file `formula_formulation.m` solves Problem 4 in the paper. Follow the instructions in that file to prove local/global optimality.

Depending on the settings you choose in each file, you may need to install [YALMIP](https://yalmip.github.io/) (tested with [this release](https://github.com/yalmip/YALMIP/releases/tag/R20200116_hotfix)), [SNOPT](https://ccom.ucsd.edu/~optimizers/) and/or [MOSEK](https://www.mosek.com/).   

> **Approval for release**: This code was approved for release by The Boeing Company in December 2020. 