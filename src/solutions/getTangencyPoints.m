% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

function tangency_points=getTangencyPoints(degree, interv)

if(degree>7 || degree<=0)
    error('not implemented yet!')
end

if(degree==1)
    tangency_points= []; 
else

sol=load(strcat('solutionTangencyPointsDeg',num2str(degree),'.mat'));
tangency_points=sol.tangencyPoints;

end


for i=1:length(tangency_points)
    tangency_points(i)=convertNumberFromABtoCD(tangency_points(i),[-1,1],interv);
end

end