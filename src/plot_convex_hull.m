function volume=plot_convex_hull(pol_x,pol_y,pol_z,A,color)
    cx=pol_x;
    cy=pol_y;
    cz=pol_z;

    vx=inv(A')*cx;
    vy=inv(A')*cy;
    vz=inv(A')*cz;

    v1=[vx(1) vy(1) vz(1)]';
    v2=[vx(2) vy(2) vz(2)]';  
    v3=[vx(3) vy(3) vz(3)]';
    v4=[vx(4) vy(4) vz(4)]';  

       
    plot3(v1(1),v1(2),v1(3),'-o','Color',color,'MarkerSize',10)
    plot3(v2(1),v2(2),v2(3),'-o','Color',color,'MarkerSize',10)
    plot3(v3(1),v3(2),v3(3),'-o','Color',color,'MarkerSize',10)
    plot3(v4(1),v4(2),v4(3),'-o','Color',color,'MarkerSize',10)
    
    [k1,volume] = convhull(vx,vy,vz);

    trisurf(k1,vx,vy,vz,'FaceColor',color)
   
    xlabel('x')
    ylabel('y')
    zlabel('z')
    alpha 0.2

end
