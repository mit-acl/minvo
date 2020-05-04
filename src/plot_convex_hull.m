
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

       
%     plot3(v1(1),v1(2),v1(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v2(1),v2(2),v2(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v3(1),v3(2),v3(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v4(1),v4(2),v4(3),'-o','Color',color,'MarkerSize',10)
    
    [x,y,z] = sphere(50);
    radius=0.02;
    color_vertex=[.98 .45 .02];
    s1 = surf(x*radius+v1(1),y*radius+v1(2),z*radius+v1(3),'FaceColor',color_vertex,'LineStyle','none' );
    s2 = surf(x*radius+v2(1),y*radius+v2(2),z*radius+v2(3),'FaceColor',color_vertex,'LineStyle','none');
    s3 = surf(x*radius+v3(1),y*radius+v3(2),z*radius+v3(3),'FaceColor',color_vertex,'LineStyle','none');
    s4 = surf(x*radius+v4(1),y*radius+v4(2),z*radius+v4(3),'FaceColor',color_vertex,'LineStyle','none');
%     
%     
%     set(s1,'facecolor',color_vertex);
%     set(s2,'facecolor',color_vertex);
%     set(s3,'facecolor',color_vertex);
%     set(s4,'facecolor',color_vertex);
    
       alpha(s1,1.0)
    alpha(s2,1.0)
    alpha(s3,1.0)
    alpha(s4,1.0)
    
%     shading faceted
%     hold on
     axis equal
     camlight
     lighting phong
%    
    [k1,volume] = convhull(vx,vy,vz);

    s2=trisurf(k1,vx,vy,vz,'LineWidth',1,'FaceColor',color)
   
    xlabel('$x$')
    ylabel('y')
    zlabel('z')
    alpha(s2,0.1)

end
