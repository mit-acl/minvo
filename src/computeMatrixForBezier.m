function Abz=computeMatrixForBezier(deg,interval)

syms t real

Abz=[];
tmp=bernsteinMatrix(deg, t);
for i=1:length(tmp)
    Abz=[Abz; double(coeffs(tmp(i),t,'All'))];

end




syms tt real

if(interval=="m11") %[-1,1]
    tt=t/2.0+0.5;
    T=[];
    
    for i=0:(deg)
        T=[tt^i T];
    end
    
    tmp=Abz*T';

    Abz=[];

    for i=1:(deg+1)
        Abz=[Abz ;double(vpa(coeffs(tmp(i),'All')))];
    end
    

    
%     subs(Abz,t,tt);
elseif(interval=="01")%[0,1]
    %Don't do anything
else
    error("not implemented yet")
end
    
end