%% For n=3

%Everything here is for the interval [0,1]
A_be=computeMatrixForBezier(3,"01");
A_bs=computeMatrixForBSpline(3,"01");
A_mv=getSolutionA(3,"01");

Mbs2mv=A_bs*inv(A_mv);


Mbs2be=A_bs*inv(A_be);

% 1/960*[182 685 100 -7;
%       56  640  280  -16;
%       -16  280  640  56;
%       -7 100 685 182;]




%% For n=2

%Everything here is for the interval [0,1]
A_be=computeMatrixForBezier(2,"01");
A_bs=computeMatrixForBSpline(2,"01");
A_mv=getSolutionA(2,"01");

Mbs2mv=A_bs*inv(A_mv);


Mbs2be=A_bs*inv(A_be);

%% For n=1

%Everything here is for the interval [0,1]
A_be=computeMatrixForBezier(1,"01");
A_bs=computeMatrixForBSpline(1,"01");
A_mv=getSolutionA(1,"01");

Mbs2mv=A_bs*inv(A_mv);


Mbs2be=A_bs*inv(A_be);