function  P=generateRandPol(deg,interv)

a = -1; b = 1;
tt=linspace(min(interv),max(interv),deg+1);
P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];

end