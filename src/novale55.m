clc; clear; close all;
t = sdpvar(1);


[p0,a0,v0] = polynomial(t,3);
[p1,a1,v1] = polynomial(t,3);
[p2,a2,v2] = polynomial(t,3);
[p3,a3,v3] = polynomial(t,3);

A=[a0';a1';a2';a3'];

Model=[];

order=15;

[s0,c0] = polynomial(t,order);
Model = [Model, sos(s0), sos(p0 - s0*(1-t^2))];

[s1,c1] = polynomial(t,order);
Model = [Model, sos(s1), sos(p1 - s1*(1-t^2))];

[s2,c2] = polynomial(t,order);
Model = [Model, sos(s2), sos(p2 - s2*(1-t^2))];

[s3,c3] = polynomial(t,order);
Model = [Model, sos(s3), sos(p3 - s3*(1-t^2))];


Model = [Model, sos(s3), a0+a1+a2+a3==[0 0 0 1]'];


settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'showprogress',1,'solver','kktqp','verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);

sol=solvesos(Model, -det(A,'polynomial'),settings,[a0;a1;a2;a3;c0;c1;c2;c3]);

% solvesos(Model, sum(a0),[],[a0;a1;a2;a3;c0;c1;c2;c3]);