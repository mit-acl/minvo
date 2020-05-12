% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

function Abz=computeMatrixForBezier(deg,interval)

syms t real

Abz=[];
tmp=bernsteinMatrix(deg, t);
for i=1:length(tmp)
    Abz=[Abz; double(coeffs(tmp(i),t,'All'))];

end

syms tt real

if(interval=="m11") %[-1,1]

   Abz= convertAFrom00toM11(Abz);

elseif(interval=="01")%[0,1]
    %Don't do anything
else
    error("not implemented yet")
end
    
end