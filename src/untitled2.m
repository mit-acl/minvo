syms t
interv=[-1,1]
x=t*t;
y=t;
z=sym(1);
figure; hold on;
fplot3(x,y,z,interv)
xlabel('x'); ylabel('y'); zlabel('z')

for t=-1:0.01:1

    x=[0 t*t];
    y=[0 t];
    z=[0 1];
    plot3(x,y,z)
end

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

axis equal

% A=0.1*[4 3   6; 
%       -3 2   -4;
%       -1 -5  -1];

  A=sol2.A;
  
for t=-1:0.01:1

    x=[0 t*t];
    y=[0 t];
    z=[0 1];
    transformed=A*[t*t,t,1]';
    
    x=[0 transformed(1)];
    y=[0 transformed(2)];
    z=[0 transformed(3)];
    
    plot3(x,y,z,'r')
end



    transformed=A*[1,0,0]';
    
    x=[0 transformed(1)];
    y=[0 transformed(2)];
    z=[0 transformed(3)];
    
    plot3(x,y,z,'b','LineWidth',3)


    transformed=A*[0,1,0]';
    
    x=[0 transformed(1)];
    y=[0 transformed(2)];
    z=[0 transformed(3)];
    
    plot3(x,y,z,'b','LineWidth',3)


    transformed=A*[0,0,1]';
    
    x=[0 transformed(1)];
    y=[0 transformed(2)];
    z=[0 transformed(3)];
    
    plot3(x,y,z,'b','LineWidth',3)