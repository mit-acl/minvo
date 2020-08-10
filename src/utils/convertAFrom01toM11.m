% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


function A_converted=convertAFrom01toM11(A)
    
    syms t
    
    tt=t/2.0+0.5;
    T=[];
    deg=size(A,1)-1;
    for i=0:(deg)
        T=[tt^i T];
    end
    
    tmp=A*T';

    A_converted=[];

    for i=1:(deg+1)
        A_converted=[A_converted ;double(vpa(coeffs(tmp(i),'All')))];
    end

end