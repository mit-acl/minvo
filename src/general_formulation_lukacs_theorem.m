% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

clear; close all; %clc;
set(0,'DefaultFigureWindowStyle','docked') %'normal'
addpath(genpath('./utils'));
addpath(genpath('./solutions'));

%delete(gcp('nocreate')); %This deletes the parallel pool

n=3;
deg_is_even = (rem(n, 2) == 0);

interv=[-1,1];

constraints=[];

use_yalmip=false; %if false, it will attempt to find a local minima

   
if(use_yalmip)
    G=[]; H=[];   
    for i=1:(n+1)
       G=[G ; sdpvar(1,floor(n/2)+1)];
       H=[H ; sdpvar(1,floor((n-1)/2)+1)];
    end
    t=sdpvar(1,1); %This is not optimized, but needed to handle symbolic expressions
else
    G=sym('G_%d%d',[n+1,floor(n/2)+1],'real');
    H=sym('H_%d%d',[n+1,floor((n-1)/2)+1],'real');
    syms t real  %This is not optimized, but needed to handle symbolic expressions
end

A=[];
lambdas=[];
for i=1:(n+1)
    
    coeff_G=G(i,:);
    coeff_H=H(i,:);
    
    TC=[];
    for i=0:(length(coeff_G)-1)
        TC=[t^i;TC]; 
    end

    TD=[];
    for i=0:(length(coeff_H)-1)
        TD=[t^i;TD]; 
    end
    
    g=(coeff_G*TC);
    h=(coeff_H*TD);
    % See now Markov–Lukács theorem   https://www.seas.ucla.edu/~vandenbe/publications/nnp.pdf  Page 952
    
    a=min(interv);
    b=max(interv);
    if(deg_is_even==1) 
        lambdai=g^2 + (t-a)*(b-t)*h^2; 
    else
        lambdai=(t-a)*g^2 + (b-t)*h^2;
    end
    
    if(use_yalmip)
        coeffs_lambdai=flip(coefficients(lambdai,t))'; 
    else
        coeffs_lambdai=coeffs(lambdai,t,'All');
    end
    A=[A; coeffs_lambdai];
    
end
sum_col_A=sum(A);
constraints=[constraints, sum_col_A==[zeros(1,size(A,1)-1) 1]   ];

%I want to maximize the absolute value of the determinant of A
disp('Computing determinant');
obj=-computeDet(A(1:end-1,1:end-1)); %This is needed when using fmincon (det(A,'polynomial') doesnt work) 
%obj=-det(A,'polynomial');

if(use_yalmip)
    general_settings=sdpsettings('savesolveroutput',0,'savesolverinput',1,'showprogress',2,'verbose',2,'debug',1);
    settings=sdpsettings(general_settings,'usex0',1,'solver','bmibnb','bmibnb.maxiter',5e20000,'bmibnb.maxtime',5e20000);
    settings.bmibnb.uppersolver='fmincon';
    settings.bmibnb.lowersolver='fmincon';
    settings.bmibnb.lpsolver='linprog';
    settings.bmibnb.lpreduce=1;
    settings.bmibnb.absgaptol=1e-15;
    tolerance=1e-11;
    settings.fmincon.TolFunValue=tolerance;
    settings.fmincon.TolFun=tolerance;
    settings.fmincon.TolConSQP=tolerance;
    settings.fmincon.TolCon=tolerance;
    settings.fmincon.TolPCG=tolerance;
    settings.fmincon.TolProjCG=tolerance;
    result=optimize(constraints,obj,settings)
else

    GH=[G H];

    lb=[];%-ones(size(GH));
    ub=[];%ones(size(GH));
    x0 = rand(size(GH));

    opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',10000,'StepTolerance',1e-9);
    problem = createOptimProblem('fmincon','objective', ...
        @(GH_x) getObj(GH_x,GH, obj),'x0',x0,'lb',lb,'ub',ub,...
        'nonlcon',@(GH_x) getConstraints(GH_x,GH, sum_col_A),'options',opts);

    gs = GlobalSearch('Display','iter');
    ms = MultiStart('Display','iter','UseParallel',true);

    disp('Running, it usually takes some time until the parpool starts');
    [xgs,fval,exitflag,output,solsgs] = run(ms,problem,2); %8000
    obj_values=[solsgs.Fval]; %In increasing order
    exit_flags=[solsgs.Exitflag];

    if(exit_flags(1)==1)
        disp("LOCAL MINIMA FOUND")
    end

end

%%
A_value=vpa(subs(A,GH,xgs));
sum(A_value)
det(A_value)

fplot(A_value*getT(n,t),[-1,1])

% getConstraints(GH_x, GH, sum_col_A)

function result=getObj(GH_x, GH, obj)
     result=double(subs(obj,GH,GH_x));
end

function [c,ceq] =getConstraints(GH_x, GH, sum_col_A)
    c=[];  %No inequality constraints
    sum_col_A_tmp=double(subs(sum_col_A,GH,GH_x));
    
    ceq=double(sum_col_A_tmp-[zeros(1,numel(sum_col_A_tmp)-1) 1]); %A'1=e
end

%%

% 
% clear t
% disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'
% 
% % check(constraints)
% %%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%To prove local optimality of n=4, uncomment this line:
% % settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'showprogress',1,'verbose',2,'debug',1);
% settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1); %snopt
% % [W_tmp, V_tmp]=findWVgivenA(getA_MV(deg,[-1,1]));
% % assign(W,W_tmp)
% % assign(V,V_tmp)
% % assign(A,getA_MV(deg,interv))
% % check(constraints)
% 
% assign(C,1/sqrt(2*(deg+1))*ones(size(C)))
% assign(D,1/sqrt(2*(deg+1))*ones(size(D)))
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%To prove global optimality of n=1,2,3, uncomment these lines:
% % constraints=[constraints obj<=-0.2]; %Needed to prevent the solver from getting stuck at A=zeros()
% % general_settings=sdpsettings('savesolveroutput',0,'savesolverinput',1,'showprogress',2,'verbose',2,'debug',1);
% % settings=sdpsettings(general_settings,'usex0',1,'solver','bmibnb','bmibnb.maxiter',5e20000,'bmibnb.maxtime',5e20000);
% % settings.bmibnb.uppersolver='fmincon';
% % settings.bmibnb.lowersolver='snopt';
% % settings.bmibnb.lpsolver='linprog';
% % settings.bmibnb.lpreduce=1;
% % settings.bmibnb.absgaptol=1e-15;
% % tolerance=1e-11;
% % settings.fmincon.TolFunValue=tolerance;
% % settings.fmincon.TolFun=tolerance;
% % settings.fmincon.TolConSQP=tolerance;
% % settings.fmincon.TolCon=tolerance;
% % settings.fmincon.TolPCG=tolerance;
% % settings.fmincon.TolProjCG=tolerance;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Other possible tests and solvers you may want to use
% % settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% % settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1);
% % %,'ipopt.tol',1e-10  %'penlab.max_outer_iter',100000
% % settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
% % result=solvemoment(constraints,obj,[],4);
% % check(constraints)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% result=optimize(constraints,obj,settings)
% 
% 
% A_minvo=value(A);
% 
% disp("abs(|A_minvo|)=")
% abs(det(A_minvo))
% 
% 
% syms t real
% T=[];
% for i=0:(deg)
%    T=[t^i ;T];
% end
%     

% fplot(A_minvo*T, [-1,1])
