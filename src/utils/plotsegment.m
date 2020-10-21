function plotsegment(a,b,color,linewidth)
   AB=[a b];
   plot3(AB(1,:),AB(2,:),AB(3,:),color,'LineWidth',linewidth)
end