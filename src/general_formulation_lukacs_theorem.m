% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'
%Useful to plot the result: http://nurbscalculator.in/
%READ THIS: https://yalmip.github.io/example/nonconvexquadraticprogramming/

deg=5;
deg_is_even = (rem(deg, 2) == 0);

% sdpvar t %This is not optimized, but needed to handle symbolic expressions
t=sdpvar(1,1);

constraints=[];

all_coeff_C=[]; all_coeff_D=[];
for i=1:(deg+1)
   all_coeff_C=[all_coeff_C ; sdpvar(1,floor(deg/2)+1)];
   all_coeff_D=[all_coeff_D ; sdpvar(1,floor((deg-1)/2)+1)];
end

A=[];
lambdas=[];
for i=1:(deg+1)
   
    
    coeff_C=all_coeff_C(i,:);
    coeff_D=all_coeff_D(i,:);
    
    TC=[];
    for i=0:(length(coeff_C)-1)
        TC=[t^i;TC]; 
    end

    TD=[];
    for i=0:(length(coeff_D)-1)
        TD=[t^i;TD]; 
    end
    
    g=(coeff_C*TC);
    h=(coeff_D*TD);
    % See now Markov–Lukács theorem   https://www.seas.ucla.edu/~vandenbe/publications/nnp.pdf  Page 952
    if(deg_is_even==1) 
        lambdai=(g)^2 + (t+1)*(1-t)*(h)^2; 
    else
        lambdai=(t+1)*(g)^2 + (1-t)*(h)^2;
    end
    coeffs_lambdai=flip(coefficients(lambdai,t))';
    A=[A; coeffs_lambdai]; 
    
end

sum_colums_A=sum(A);
constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
constraints=[constraints, sum_colums_A(end)==1];

%I want to maximize the absolute value of the determinant of A
obj=-computeDet(A(1:end-1,1:end-1)); %This is needed when using fmincon (det(A,'polynomial') doesnt work) 
%obj=-det(A,'polynomial');
%%


clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'

% check(constraints)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%To prove local optimality of n=4, uncomment this line:
% settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'showprogress',1,'verbose',2,'debug',1);
settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1);
% [W_tmp, V_tmp]=findWVgivenA(getA_MV(deg,[-1,1]));
% assign(W,W_tmp)
% assign(V,V_tmp)
% assign(A,getA_MV(deg,[-1,1]))
% check(constraints)

assign(all_coeff_C,1/sqrt(2*(deg+1))*ones(size(all_coeff_C)))
assign(all_coeff_D,1/sqrt(2*(deg+1))*ones(size(all_coeff_D)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%To prove global optimality of n=1,2,3, uncomment these lines:
% constraints=[constraints obj<=-0.2]; %Needed to prevent the solver from getting stuck at A=zeros()
% general_settings=sdpsettings('savesolveroutput',0,'savesolverinput',1,'showprogress',2,'verbose',2,'debug',1);
% settings=sdpsettings(general_settings,'usex0',1,'solver','bmibnb','bmibnb.maxiter',5e20000,'bmibnb.maxtime',5e20000);
% settings.bmibnb.uppersolver='fmincon';
% settings.bmibnb.lowersolver='snopt';
% settings.bmibnb.lpsolver='linprog';
% settings.bmibnb.lpreduce=1;
% settings.bmibnb.absgaptol=1e-15;
% tolerance=1e-11;
% settings.fmincon.TolFunValue=tolerance;
% settings.fmincon.TolFun=tolerance;
% settings.fmincon.TolConSQP=tolerance;
% settings.fmincon.TolCon=tolerance;
% settings.fmincon.TolPCG=tolerance;
% settings.fmincon.TolProjCG=tolerance;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Other possible tests and solvers you may want to use
% settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1);
% %,'ipopt.tol',1e-10  %'penlab.max_outer_iter',100000
% settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
% result=solvemoment(constraints,obj,[],4);
% check(constraints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


result=optimize(constraints,obj,settings)


A_minvo=value(A);

disp("abs(|A_minvo|)=")
abs(det(A_minvo))


syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
    

% fplot(A_minvo*T, [-1,1])
