
%Using robust optimization, doesn't seem to work for my problem

%Read this tutorial for more details: https://yalmip.github.io/tutorial/robustoptimization/
clc; cloar; close all;

n=3;
A=sdpvar(n+1, n+1,'full');
sdpvar t

constraints=[sum(A)==[zeros(1,size(A,1)-1) 1] ];
constraints = [constraints -1.0<= t <= 1.0, uncertain(t)];
constraints=[constraints A*getT(n,t)>=zeros(n+1,1) ];

settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);


A_cropped=A(1:end-1,1:end-1);
obj=-det(A_cropped,'polynomial');
sol = optimize(constraints,obj,settings)