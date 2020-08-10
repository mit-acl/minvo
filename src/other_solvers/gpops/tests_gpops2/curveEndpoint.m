% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


function output = curveEndpoint (input)

[B R]=generateBR(input.phase.finalstate);
A=getA(B,R);

obj=abs(double(det(A)));

output.objective = -obj; %I want to maximize |A|



