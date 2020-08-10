from ncpol2sdpa import *
import numpy as np
from sympy import *

level = 3


deg=3;
num_el_A_cropped=(deg)*(deg+1);
A_cropped_var = generate_variables('A', num_el_A_cropped)
print A_cropped_var
A_cropped=np.reshape(A_cropped_var, (deg, deg+1))
columns_sum=A_cropped.sum(axis=0);
print((np.array([0,0,0,1])-columns_sum)).shape
tmp=np.reshape([0,0,0,1]-columns_sum,(1,deg+1));
A=np.concatenate((A_cropped, tmp));

print(A)


t=[-1.0, -0.77354858731953, -0.030882541230, 0.030882541230, 0.77354858731953, 1.0];

inequalities=[]
for t_i in t:
    T=np.array([t_i**3, t_i**2, t_i**1, 1])
    T=np.reshape(T, (deg+1, ));
    for j in range(deg+1):
        row_of_A=np.reshape(A[j:j+1,:], (deg+1, ));
        inequalities.append(np.inner(row_of_A,T))

print(inequalities)

#Convert to simpy matrix
A = Matrix(A)
detA=A.det()

obj=-detA;

sdp=SdpRelaxation(flatten([A_cropped_var]),verbose=1)
sdp.get_relaxation(level, objective=obj, inequalities=inequalities, chordal_extension=True)
sdp.solve()
print(sdp.primal, sdp.dual)