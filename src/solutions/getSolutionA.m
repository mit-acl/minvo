% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

function [A rootsA]=getSolutionA(degree, interval)

%Note that these guesses are in [0,1]
sol=load(strcat('solutionDeg',num2str(degree),'.mat'));
A_saved=sol.A;

sol=load(strcat('rootsLambdaiDeg',num2str(degree),'.mat'));
roots_lambdai=sol.roots_lambda_solution;

% order stuff
if(degree==1)
    co=[1 2]; %correct order
elseif(degree==2)
    co=[3 2 1]; %correct order
elseif(degree==3)
    co=[1 4 2 3]; %correct order    
elseif(degree==4)
    co=[5 4 3 1 2]; %correct order
elseif(degree==5)
    co=[3 4 2 5 1 6]; %correct order
elseif(degree==6)
    co=[2 5 3 4 7 1 6]; %correct order
elseif(degree==7)
    co=[3 6 1 8 4 5 2 7]; %correct order
end




A=[];
rootsA=[];
i=1;
for j=co
    A=[A; A_saved(j,:)];
    rootsA{i}=roots_lambdai{j};
    i=i+1;
end


%Sort the roots in increasing order
for j=1:length(rootsA)
    rootsA{j}=reshape(rootsA{j},1,[]); %Make sure its a row vector;
    rootsA{j}=sort(rootsA{j});
end

if(interval=="01")
    A=convertAFromM11to01(A);
    
    for j=1:length(rootsA)
        rootsA{j}=convertFromM11to01(rootsA{j});
    end
    
end



end


function x_converted=convertFromM11to01(x)
   x_converted=(x+1)/2.0; %This applies to all the vector
end