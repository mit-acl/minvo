function plotCurve(P,interv)
    deg=size(P,2)-1; %degree
    syms t real
    T=getT(deg,t);  fplot(P(1,:)*T,P(2,:)*T,interv,'r','LineWidth',2);
end