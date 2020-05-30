
function Abs=computeMatrixForBSpline(deg,interval)



if(deg==1)
    Abs=[-1  1 ;
         1  0 ];
elseif(deg==2)
    Abs=(1/2)*[1  -2   1 ;
            -2   2  1 ;
             1   0   0  ];
elseif(deg==3)
    Abs=(1/6)*[-1 3 -3 1;
            3 -6 0 4;
            -3 3 3 1;
            1 0 0 0 ];
else
    error("Not implemented yet")
end


if(interval=="m11") %[-1,1]

   Abs= convertAFrom01toM11(Abs);

elseif(interval=="01")%[0,1]
    %Don't do anything
else
    error("not implemented yet")
end
    
end