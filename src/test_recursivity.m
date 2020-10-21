%%This file tries to check numerically if there exists a recursive formula of the polynomials
%%of the MINVO basis. It does so by imposing conditions between the
%%polynomials of different n, and then solving by least squares

close all; clear; clc;

addpath(genpath('./utils'));
addpath(genpath('./solutions'));

global interv
interv=[0,1];

t=sym('t',[1,1], 'real')
getLambda(2,1,t)

n_initial=2;
n_final=2;
q=9;

a=sym('a',[1,q+1],'real');
b=sym('b',[1,q+1],'real');
c=sym('c',[1,q+1],'real');

eqs=[];


for n=n_initial:n_final
   
   for i=1:(n+1)    
       
       lambda_np1_i= getLambda(n+1,i,t);
       lambda_n_i= getLambda(n,i,t);
       lambda_n_im1= getLambda(n,i-1,t);
       
       blend1=a*getT(q,t);
       blend2=b*getT(q,t);

       should_be_zero=lambda_np1_i-(blend1*lambda_n_i+blend2*lambda_n_im1 );

%        should_be_zero=lambda_np1_i-(blend1*lambda_n_i+blend2*lambda_n_im1);
       
       coeff1=coeffs(should_be_zero,t,'All');
       eqs=[eqs coeff1==zeros(size(coeff1))];
   end
    
end


[C,d]=equationsToMatrix(eqs);
C=double(C); d=double(d); %% Cx=d

[x,resnorm] =  lsqlin(C,d,[],[])

if(resnorm<0.1) %If the residual is small enogh 
    disp("Recursive formula Found")
else
    disp('NOT found')
end
    


function result=getLambda(n,i,t)
  global interv
  
  
%   if(i<=0)
%       i
%       i=i+n+1
%   end
%   
%   if(i>(n+1))
%       i
%       i=i-(n+1)
%   end
%   disp('===')
%   A=getA_MV(n,interv);
%    result= A(i,:)*getT(n,t);
     
  
  if(i<1 || i>(n+1))
      result=0*t;
  else
      A=getA_MV(n,interv);
      result= A(i,:)*getT(n,t);
  end
  


end