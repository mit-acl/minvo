# MINVO Basis: Finding Simplexes with Minimum Volume Enclosing Polynomial Curves #

![](./imgs/minvo3d.png) 

## Citation

When using MINVO, please cite [this paper](https://www.google.com/):

```bibtex
@inproceedings{tordesillas2019faster,
  title={{MINVO} Basis: Finding Simplexes with Minimum Volume Enclosing Polynomial Curves},
  author={Tordesillas, Jesus and How, Jonathan P},
  journal={arXiv preprint},
  year={2020}
}
```

## Instructions to use the MINVO basis

### For a particular curve (Problem 1 of the paper)
Given a matrix P that contains in each row

### For a particular simplex  (Problem 2 of the paper)
Given the matrix <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{V}"> whose columns contain the vertexes of the <img src="https://render.githubusercontent.com/render/math?math=n">-simplex, you can obtain the <img src="https://render.githubusercontent.com/render/math?math=n">-th order polynomial curve that is enclosed in that simplex that has largest convex hull by running this:

```matlab
addpath(genpath('./solutions'));
P=V*getSolutionA(n,"m11");
```
Each row <img src="https://render.githubusercontent.com/render/math?math=i"> of <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{P}"> will contain the coefficients (in decreasing order) of the polynomial of the <img src="https://render.githubusercontent.com/render/math?math=i">-th coordinate of the curve <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{p}(t) t \in [-1,1]">.


## Instructions to derive the MINVO basis

### To generate all the figures of the paper
You can simply run `plot_solution.m` (you don't need to install any external solver to do this).

### To run the optimization
To 
