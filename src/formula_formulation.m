% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

close all; clc; clear;
%delete(gcp('nocreate')); %Delete the parallel pool
global use_yalmip

%%%%%%%% Instructions:
% Set use_yalmip=true to prove global optimality for n=1,2,3 (using moment relaxations)
% Set use_yalmip=false to prove local optimality for n=1,2,3,4,5,6,7 (using fmincon from many different initial guesses)
use_yalmip=false;
n=3;

%
if (use_yalmip==true)
    sdpvar t
else
    syms t
end

deg=n;
deg_is_even = (rem(deg, 2) == 0);

W=[]; B=[];R=[];
if(deg_is_even==0) %Deg is odd
    
    jj=1;
    for i=1:2:deg
       
        B=appendElement(B,strcat('B',num2str(i)));  Bi=B(end);
        pol=-Bi*(t-1); roots_lambda{jj}=[1.0];
        for j=1:(deg-1)/2                                
            R=appendElement(R,strcat('R',num2str(i),num2str(j)));  Rij=R(end);          
            pol=pol*((t-Rij)^2);
        end
        W=[W;pol];
        jj=jj+1;
    end
    
    W=[W;substitute(W,t,-t)]; %Insert the other half
    
else %Deg is even

    n_2_is_even = (rem(n/2, 2) == 0);
    if(n_2_is_even)
        Ia=2:2:(n/2);
        Ib=1:2:(n/2); %n/2+1
    else
        Ia=2:2:(n/2); %n/2+1
        Ib=1:2:(n/2); 
    end
    

    jj=1;
    for index=1:size(Ia,2)
        i=Ia(index);
        B=appendElement(B,strcat('B',num2str(i)));  Bi=B(end);
        pol=-Bi*(t+1)*(t-1); roots_lambda{jj}=[-1.0, 1.0];
        for j=1:(n-2)/2
            R=appendElement(R,strcat('R',num2str(i),num2str(j)));  Rij=R(end);
            pol=pol*((t-Rij)^2);
        end
        W=[W;pol];
        jj=jj+1; 
    end  
 
   for index=1:size(Ib,2)
       i=Ib(index);
       B=appendElement(B,strcat('B',num2str(i)));  Bi=B(end);
       pol=Bi;   roots_lambda{jj}=[];
       for j=1:(n/2)
            R=appendElement(R,strcat('R',num2str(i),num2str(j)));  Rij=R(end);
            pol=pol*((t-Rij)^2);
       end
       W=[W;pol];
       jj=jj+1;
   end  
   
   
   %Force the lambda in the middle to have symmetrical roots
   i=n/2+1;  jj=i;
   if(n_2_is_even)
       B=appendElement(B,strcat('B',num2str(i)));  Bi=B(end);
       pol=B(i);   roots_lambda{jj}=[];
       for j=1:(n/4)
            R=appendElement(R,strcat('R',num2str(i),num2str(j)));  Rij=R(end);
            pol=pol*((t-Rij)^2)*((t+Rij)^2);
       end
       W=[W;pol];
   else
       B=appendElement(B,strcat('B',num2str(i)));  Bi=B(end);
       pol=-B(i)*(t+1)*(t-1); roots_lambda{jj}=[1.0]; %Adding only 1.0 because then I'll later make the symmetric. 
       for j=1:((n-2)/4)
            R=appendElement(R,strcat('R',num2str(i),num2str(j)));  Rij=R(end);
            pol=pol*((t-Rij)^2)*((t+Rij)^2);
       end
       W=[W;pol];
   end
   
   W=[W;substitute(W(1:end-1,:),t,-t)]; %Force symmetry in the first block except the last row
   
end
%%
%Get the coefficients
coeffic=getCoeffInDecOrder(sum(W), t, deg);

%% Solve using yalmip
if(use_yalmip==true)
        
    disp("Creating the A matrix")
    %Create the A matrix:
    A=[];
    for i=1:length(W)
        tmp=getCoeffInDecOrder(W(i),t,deg);
        A=[A; tmp];
    end
    
    clear t
    constraints=[];

    % NOTE: 
    % constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
    % constraints=[constraints, sum_colums_A(end)==1];
    % don't work when sum_colums_A has one value that is '0'

    sum_colums_A=sum(A);

    for i=1:(size(sum_colums_A,2)-1)
        if( double(sum_colums_A(i))==0)
            %nothing
        else
            constraints=[constraints, sum_colums_A(i)==0.0   ];
        end
    end

    constraints=[constraints, sum_colums_A(end)==1];
    constraints=[constraints, B>=zeros(size(B))];

    obj=-det(A,'polynomial');% computeDet(A)

    disp("Starting optimization")
    [sol,x,momentdata]=solvemoment(constraints,obj,[], 5); 
   
%%%Other possible settings you may want to explore
%     settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
%     assign(A,A_guess);
%     usual_settings=sdpsettings('savesolveroutput',1,'savesolverinput',1,'showprogress',1,'verbose',2,'debug',1);
%     settings=sdpsettings('usex0',0,'solver','moment',usual_settings);
%     settings=sdpsettings(usual_settings,'solver','sparsepop','sparsepop.relaxOrder',4);
%     settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
%     settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
%     settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','ipopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);


else
%% Solve using fmincon (GlobalSearch)

    %%%%%%%%%%%%%%%% This section simply stores in a struct the roots of
    % each lambda_i, to be able to recover them later
    Atmp=[];
    for i=1:length(W)
            tmp=getCoeffInDecOrder(W(i),t,deg);
            Atmp=[Atmp; tmp];
    end    
    
    for jj=1:(ceil(size(Atmp,1)/2))
        tmp=symvar(Atmp(jj,:));
        tmp=subs(tmp,B,zeros(size(B))); %set all the values of B to zero
        roots_lambda_jj=nonzeros(tmp)'; %get all the variables Rij
        roots_lambda{jj}=[roots_lambda{jj} roots_lambda_jj];
    end
    
    if(deg_is_even)
        tmp2=ceil(size(Atmp,1)/2);
        roots_lambda{tmp2}=[roots_lambda{tmp2} -roots_lambda{tmp2}]; %symmetrical roots for this one
    end
    
    %%%%%%%%%%%%%%%
    %Solve for the coefficients:
    disp("Solving linear system")
    solution=solve(coeffic==[zeros(1,deg) 1],B); %This is a linear system

    disp("struct2array")
    if(deg~=1)
        B_solved=(struct2array(solution)');
        B_solved=B_solved(:);
    else
        B_solved=solution;
    end

    disp("Substituting")
    W=subs(W,B,B_solved);

    disp("Creating the A matrix")
    %Create the A matrix:
    A=[];
    for i=1:length(W)
            tmp=getCoeffInDecOrder(W(i),t,deg);
            sdisplay(tmp)
            A=[A; tmp];
    end

    disp("Computing the determinant")
    %Compute the determinant
    if(deg_is_even==0)
        determ=computeDetSmartly(A);
    else
        determ=det(A);
    end

    lb=-ones(size(R));
    ub=ones(size(R));
    x0 = rand(size(R));

    opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',10000,'StepTolerance',1e-9);
    problem = createOptimProblem('fmincon','objective', ...
        @(R_x) getObj(R_x,R, determ),'x0',x0,'lb',lb,'ub',ub,...
        'nonlcon',@(R_x) getConstraints(R_x,R, B_solved),'options',opts);

    % Construct a GlobalSearch object
    gs = GlobalSearch('Display','iter');
    % Construct a MultiStart object based on our GlobalSearch attributes
    ms = MultiStart('Display','iter','UseParallel',true);

    disp('Running, it usually takes some time until the parpool starts');
    [xgs,fval,exitflag,output,solsgs] = run(ms,problem,400); %8000

    
    obj_values=[solsgs.Fval]; %In increasing order
    exit_flags=[solsgs.Exitflag];
    
    if(exit_flags(1)==1)
        disp("LOCAL MINIMA FOUND")
    end

    %%Recover solution
    B_solution=vpa(subs(B_solved,R,xgs));
    A_solution=double(vpa(subs(A,R,xgs)));
    R_solution=double(vpa(subs(R,R,xgs)));
    
    tangencyPoints=sort([R_solution; -R_solution]); %Also the symmetric roots
    
    for i=1:size(roots_lambda,2)
        roots_lambda_solution{i}=double(subs(roots_lambda{i}, R, xgs));
    end
    
    s=size(roots_lambda,2);
    for i=1:s
        if(i==s && deg_is_even)
            continue %this lambda doesn't have a symmetric one
        end
        roots_lambda_solution{i+s}=-roots_lambda_solution{i}; %Symmetry for the rest of lambda_i
    end
    
    det(A_solution)

    A=A_solution;

    syms t real
    T=[];
    for i=0:(deg)
       T=[t^i ;T];
    end
    interv=[-1,1];
    figure;
    fplot(A_solution*T,interv)
    
%     save(['solutionDeg' num2str(deg) '.mat'],'A');%,'rootsA'
%     save(['sols_formula/solutionTangencyPointsDeg' num2str(deg) '.mat'],'tangencyPoints');
%     save(['sols_formula/rootsLambdaiDeg' num2str(deg) '.mat'],'roots_lambda_solution');

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS


%Get the coefficients in decreasing order
function coeffic=getCoeffInDecOrder(expression, variable, deg)
    global use_yalmip
    
    if (use_yalmip==true)
        %NOTE: IF YOU USE
        %coeffic=flip(coefficients(expression,variable)'), it returns [1 4]
        %for the polynomial t*t +4 (i.e. it doesn't include a 0).
        %That's why the code below is needed. 
         [c,v] = coefficients(expression,variable);
%          sdisplay(c)
         coeffic=[];%sdpvar(1,deg+1);
         for i=0:deg
             found=false;
             for j=1:size(v,1)
                if(double(v(j)-variable^i)==0)
                  coeffic=[coeffic c(j)];
                  found=true;
                 break %leave inner-most loop
                end
             end
             
             if(found==false)
                 coeffic=[coeffic 0];
             end
         end
        
        coeffic=flip(coeffic);
%         coeffic=flip(coefficients(expression,variable)');
    else
        coeffic=coeffs(expression,variable,'All'); %When using All, no flip!! (gives directly [a  b c d]
    end
    
    coeffic=[zeros(1,deg+1-length(coeffic)) coeffic];

end




function result=getObj(R_x, R, determ)
    
%     result=det(subs(A,R,R_x));
    %global R determ
     result=subs(determ,R,R_x);
    result=-abs(double(result));
end

function [c,ceq] =getConstraints(R_x, R, B_solved)
    %global R B_solved
    c=-subs(B_solved(:),R,R_x); %  -B(i) <=0
%     c=[c ; subs(R(:),R,R_x) - 2*ones(size(R(:)))]; %R(i) <= 2
%     c=[c ; -subs(R(:),R,R_x) - 2*ones(size(R(:)))]; % R(i) >=-2 //// -2 - R(i) <=0
    c=double(c);
    ceq=[];
end

function result=computeDetSmartly(A)

%odd columns of A
oddcol=A(:,1:2:end);

%even columns 
evencol=A(:,2:2:end);

A_reorganized=[oddcol,evencol];

nrow=size(A,1);
ncol=size(A,2);

AA=A_reorganized(1:nrow/2,1:ncol/2);

BB=A_reorganized(1:nrow/2,ncol/2+1:end);

result=det(2*AA)*det(BB);

end


function result=substitute(expression, var_old, var_new)
    global use_yalmip
    if (use_yalmip==true)
        expression=replace(expression,var_old,var_new);
    else
       expression=subs(expression,var_old,var_new);
    end
    result=expression;
end

function Q_result=appendElement(Q,name)
    global use_yalmip
    if (use_yalmip==true)
        Q=[Q; sdpvar(1,1)]; 
    else
       Q=[Q; sym(name,'real')];
    end
            
    Q_result=Q;
end
