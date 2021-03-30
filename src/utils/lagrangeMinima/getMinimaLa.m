function minim=getMinimaLa(i,n,interv)

    dist=linspace(min(interv),max(interv),n+1);
    P=double(lagrangePoly(dist));

    p=P(i,:);
    dp=polyder(p);
    rdp=roots(dp);

    ddp=polyder(dp);
    ddp_eval=polyval(ddp,rdp);

    minim=rdp(ddp_eval>0);

%     minim=[min(interv); minim; max(interv)];
   
end