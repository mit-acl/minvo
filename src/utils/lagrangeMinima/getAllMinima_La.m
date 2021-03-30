%Get the minima of the Lagrange basis functions in interv
%min(interv) and max(interv) are also added to the returned value
function [all_minima, all_indexes]=getAllMinima_La(n, interv)

    dist=linspace(min(interv),max(interv),n+1);
    P=double(lagrangePoly(dist));

    all_minima=[];
    all_indexes=[];
    for i=1:size(P,2)
        p=P(i,:);
        dp=polyder(p);
        rdp=roots(dp);

        ddp=polyder(dp);
        ddp_eval=polyval(ddp,rdp);

        minim=rdp(ddp_eval>0);

        all_minima=[all_minima minim']
        all_indexes=[all_indexes i*ones(1,numel(minim))]

    end

%      all_minima=[all_minima min(interv) max(interv)];
    [all_minima, order]=sort(all_minima);
    all_indexes = all_indexes(:,order);
end