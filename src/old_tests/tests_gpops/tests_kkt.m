close all; clc; clear;

Q=[eye(3) zeros(3,1)];

R = sym('R', [4 4]);

tmp= sym('V', [4,3]);
V=tmp';

lambda=sym('L',[4,1]);

Q*R

s=[0 0 0 1]';

R'*s;

V*lambda;

lambda'*ones(4,1);

R'*s;

alp_a=sym('alp_a',[3,1]);
alp_c=sym('alp_c',[4,1]);

Q'*alp_a;

s*alp_c'*ones(4,1);

s*alp_c'*R';


lambda'*R'*s;

Q'*alp_a*lambda'*R';


% En el ejemplo, el optimo es:
syms t;
lambda1= 3.4417150113470782457625318784267*t^3 - 3.3354689765319474048510528518818*t^2 + 0.80812439015635884054233883944107*t + 0.00000099938873458182832544278821618189;
lambda2=-3.4415628781064451224835920584155*t^3 + 6.9894219682793146120047822478227*t^2 - 4.4622188271595399911007007176522*t + 0.91436017358236654217762406915426;
lambda3= -6.6793936856534701362875239283312*t^3 + 8.1920134593848583648423300473951*t^2 - 1.598257232203776689871688176936*t + 0.085637805008032086284686101862462;
lambda4= 6.6792415524128374570977939583827*t^3 - 11.845966451132227348352898843586*t^2 + 5.2523516692069565081624205049593*t + 0.0000010220208668471246456828805532213;

lambda=[lambda1;lambda2;lambda3;lambda4];

lambda=vpa(lambda);

pol_x=[0.2 0.3 2 1]';%[a b c d]
pol_y=[-0.3 +3 -5 6]';%[a b c d]
pol_z=[1 -0.1 -1 -4]';%[a b c d]

T=[t*t*t; t*t; t; 1];

m=[pol_x'*T;  pol_y'*T;  pol_z'*T];

R=   [3.6990    0.8291    2.8243    1.3755;
     3.5905    6.2290    3.5546    4.8692;
    -4.0820   -3.9524   -4.5086   -4.2921;
     1.0000    1.0000    1.0000    1.0000];


vpa(lambda'*R'*s)


constraint=vpa(2*det(R)*det(R)*inv((R))-Q'*alp_a*lambda'-s*alp_c');

subs(constraint, t, 0.5)

% subs(constraint*ones(4,1), t, 0.5)
