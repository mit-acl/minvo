% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


function [B,R]= generateBR(parameters)

    global deg_pol 

    B=parameters(1,1:((deg_pol+1)/2));
    R=parameters(1,(deg_pol+1)/2+1:end);
    
    R=reshape(R, (deg_pol+1)/2, (deg_pol-1)/2);
    
    
end