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

deg=3;
deg_is_even = (rem(deg, 2) == 0);

if(deg_is_even==1)
    d=deg/2;
    size_Wi=d+1;
    size_Vi=d;
else
    d=(deg-1)/2;
    size_Wi=d+1;
    size_Vi=d+1;
end


W=[]; V=[];
for i=1:(deg+1)
   W=[W sdpvar(size_Wi,size_Wi)];
   V=[V sdpvar(size_Vi,size_Vi)];
end


sdpvar t %This is not optimized, but needed to handle symbolic expressions

constraints=[];

A=[];
lambdas=[];
for i=1:(deg+1)
    Wi=W(:,(i-1)*size_Wi+1:i*size_Wi);
    Vi=V(:,(i-1)*size_Vi+1:i*size_Vi);

    %Wi and Vi are psd matrices <=> All ppal minors are >=0
    constraints=[constraints, getPsdConstraints(Wi), getPsdConstraints(Vi)];   
%   constraints=[constraints, Wi>=0, Vi>=0]; %Works when using solvemoment sdp
    
    Tvi=[];
    for i=0:(size(Vi,1)-1)
        Tvi=[Tvi; t^i]; 
    end

    Twi=[];
    for i=0:(size(Wi,1)-1)
        Twi=[Twi; t^i]; 
    end

    if(deg_is_even==1)
        lambdai=Twi'*Wi*Twi + (t+1)*(1-t)*Tvi'*Vi*Tvi;
    else
        lambdai=(t+1)*Twi'*Wi*Twi + (1-t)*Tvi'*Vi*Tvi;
    end
    coeffs_lambdai=flip(coefficients(lambdai,t))';
    A=[A; coeffs_lambdai]; 
    
end

sum_colums_A=sum(A);
constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
constraints=[constraints, sum_colums_A(end)==1];

%I want to maximize the absolute value of the determinant of A
obj=-computeDet(A); %This is needed when using fmincon (det(A,'polynomial') doesnt work) 
%obj=-det(A,'polynomial');



clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'

% check(constraints)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%To prove local optimality of n=4, uncomment this line:
%settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1);
[W_tmp, V_tmp]=findWVgivenA(getSolutionA(deg,"m11"));
assign(W,W_tmp)
assign(V,V_tmp)
%%%% assign(A,getSolutionA(deg,"m11"))
check(constraints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%To prove global optimality of n=1,2,3, uncomment these lines:
constraints=[constraints obj<=-0.2]; %Needed to prevent the solver from getting stuck at A=zeros()
general_settings=sdpsettings('savesolveroutput',0,'savesolverinput',1,'showprogress',2,'verbose',2,'debug',1);
settings=sdpsettings(general_settings,'usex0',1,'solver','bmibnb','bmibnb.maxiter',5e20000,'bmibnb.maxtime',5e20000);
settings.bmibnb.uppersolver='fmincon';
settings.bmibnb.lowersolver='fmincon';
settings.bmibnb.lpreduce=1;
settings.bmibnb.absgaptol=1e-15;
tolerance=1e-11;
settings.fmincon.TolFunValue=tolerance;
settings.fmincon.TolFun=tolerance;
settings.fmincon.TolConSQP=tolerance;
settings.fmincon.TolCon=tolerance;
settings.fmincon.TolPCG=tolerance;
settings.fmincon.TolProjCG=tolerance;
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
    

fplot(A_minvo*T, [-1,1])

function constraints=getPsdConstraints(A)
constraints=[];

%See Theorem 17 Of Lecture 5 of https://homepages.cwi.nl/~monique/eidma-seminar-parrilo/Parrilo-LectureNotes-EIDMA.pdf

sdpvar lambda;
n=size(A,1);
charac_poly=computeDet(lambda*eye(size(A,1))-A);

coeffs_charac_poly=coefficients(charac_poly,lambda); %[p0 p1 ... p_n-1 1]

size(coeffs_charac_poly)
for(i=0:(n-1))
    j=i+1;
    constraints=[constraints coeffs_charac_poly(j)*((-1)^(n-i))>=0 ];
end
% 
% if(size(A,1)==1)
%     constraints=[A(1,1)>=0];
%     return;
% elseif(size(A,1)==2)
%     constraints=[A(1,1)>=0, A(2,2)>=0, (A(1,1)*A(2,2)-A(1,2)*A(2,1))>=0]  
%     return;
% elseif(size(A,1)==3)
% 
%     tmp=A(2:3, 2:3); %delete 1
%     constraints=[constraints getPsdConstraints(tmp)] 
%     tmp=A([1, 3:end], [1, 3:end]); %delete 2
%     constraints=[constraints getPsdConstraints(tmp)]
%     tmp=A(1:2, 1:2); %delete 3
%     constraints=[constraints getPsdConstraints(tmp)] 
%     
%     constraints=[constraints computeDet(A)>=0]
%     return;
% else
%     error('NOT YET IMPLEMENTED')
% end
    
   
end