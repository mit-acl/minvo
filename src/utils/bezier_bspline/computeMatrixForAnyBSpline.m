function A=computeMatrixForAnyBSpline(deg, index_t_start, knots,interval)

M=getM(deg, index_t_start, knots,interval);

%And now change the order of the rows/columns to the convention I use
A=[];
for i=1:size(M,1) 
    A=[M(i,:)'  A];
end


A=convertCoeffMatrixFromABtoCD(A,[0,1],interval); 

end

function M=getM(deg, index_t_start, knots,interval)

%index_t_start is 1-based

%Everything here is for interval t \in [0,1]

%Following the notation from
%https://link.springer.com/article/10.1007/s003710050206 
%("General matrix representations for B-splines"
% See Theorem 1 of page 180


if(deg>0)
    Mkm1=getM(deg-1, index_t_start, knots,interval);

    
    k=deg+1;
    i=index_t_start-1;
   

    d0_vector=[];
    for jj=0:(k-2)
        j=i-k+2+jj;
        d0_vector=[d0_vector getd0(i,j,k,knots)]; 
    end
       
   d1_vector=[];
   for jj=0:(k-2)
        d1_vector=[d1_vector getd1(i,i-k+2+jj,k,knots)]; 
   end
    
    Mk=[Mkm1; zeros(1,size(Mkm1,2))]*[diag([1-d0_vector])+(k>2)*diag([d0_vector(1:end-1)],1)  [zeros(k-2,1);d0_vector(end)]]+...
      +[zeros(1,size(Mkm1,2)); Mkm1]*[diag([-d1_vector])+(k>2)*diag([d1_vector(1:end-1)],1)  [zeros(k-2,1);d1_vector(end)]]   ;
    
    M=Mk;
else
    M=1;
end




end

function result=getd0(i,j,k,knots)
    result=(knots(i+1)-knots(j+1))/(knots(j+k-1+1)-knots(j+1)); %Note that +1 is added due to the Matlab 1-indexing
end

function result=getd1(i,j,k,knots)
    result=(knots(i+1+1)-knots(i+1))/(knots(j+k-1+1)-knots(j+1));%Note that +1 is added due to the Matlab 1-indexing
end

%    Mkm1=[1 0;-1 1];
%    Mkm1=[ 0.5000    0.5000         0;
%           -1.0000    1.0000         0;
%            0.5000   -1.0000    0.5000];

%  clc; clear;
% computeMatrixForAnyBSpline_buena(4,3,[0 0 0 0 1 1 1 1],[0,1])
% 
% deg=3;
% index_t_start=4;
% knots=[0 1 2 3 4 5 6];
% interval=[0,1];