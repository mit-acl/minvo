function  P=generateRandSmoothPol(deg,interv)

a = -1; b = 1;
tt=linspace(min(interv),max(interv),deg+1);


all_xy=zeros(2,1);

amplitude=1.0;%(b-a).*rand(2,1) + a;
for i=1:(numel(tt)-1)
    xy_new=all_xy(:,end)+rand(2,1)-0.15;
    all_xy=[all_xy xy_new];
end

x=all_xy(1,:);
y=all_xy(2,:);

P= [polyfit(tt,x,deg); 
    polyfit(tt,y,deg)];


end