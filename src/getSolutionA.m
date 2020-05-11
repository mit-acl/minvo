function A=getSolutionA(degree, interval)

%Note that these guesses are in [0,1]
guess=load(strcat('sols_formula/solutionDeg',num2str(degree),'.mat'));
A=guess.A;

% order stuff
if(degree==2)
    tmp1=A(1,:);
    tmp3=A(3,:);
    A(1,:)=tmp3;
    A(3,:)=tmp1;
elseif(degree==3)
    tmp1=A(1,:);
    tmp2=A(2,:);
    tmp3=A(3,:);
    tmp4=A(4,:);
    A=[tmp2; tmp3; tmp1; tmp4];
elseif(degree==4)
    tmp1=A(1,:);
    tmp2=A(2,:);
    tmp3=A(3,:);
    tmp4=A(4,:);
    tmp5=A(5,:);
    A=[tmp3; tmp2; tmp5; tmp1; tmp4];
elseif(degree==5)
    tmp1=A(1,:);
    tmp2=A(2,:);
    tmp3=A(3,:);
    tmp4=A(4,:);
    tmp5=A(5,:);
    tmp6=A(6,:);
    A=[tmp3; tmp4; tmp2; tmp5; tmp1; tmp6];
    % sol5.A=[tmp4; tmp1; tmp5; tmp3; tmp2; tmp6];
elseif(degree==6)
    tmp1=A(1,:);
    tmp2=A(2,:);
    tmp3=A(3,:);
    tmp4=A(4,:);
    tmp5=A(5,:);
    tmp6=A(6,:);
    tmp7=A(7,:);
     A=[tmp4; tmp1; tmp6; tmp7; tmp5; tmp2; tmp3];
elseif(degree==7)
    tmp1=A(1,:);
    tmp2=A(2,:);
    tmp3=A(3,:);
    tmp4=A(4,:);
    tmp5=A(5,:);
    tmp6=A(6,:);
    tmp7=A(7,:);
    tmp8=A(8,:);
    A=[tmp1;tmp7;tmp4;tmp6;tmp2;tmp8;tmp3;tmp5;];
end



if(interval~="m11")
    error("TODO")
end

end