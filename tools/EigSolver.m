clear;clc;
syms rho u Es Y positive
syms A1 qcon c_m positive

rho = double(0.778319266566257);
u = 1.000000011686097E-8;
Es = 1.22269771196144;
Y = 1.00000000000000;
A1 = 1.00000000000000;
qcon = 0.218110130028006;

m1 = rho;
m2 = rho*u;
m3 = rho*Es;
m4 = rho*Y;
c_m = (u^2 - 2*A1*Es + 2*A1*Y*qcon)^(1/2);

EVR = [
[ (2*A1*m1^2)/m2, -(2*A1*m1^2*qcon)/m2,                                                                                                                                   m1/m4,                                                                                                                                   m1/m4];
[     (2*A1*m1) ,     -(2*A1*m1*qcon) ,                                                                   (3*m2 + 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2))/(3*m4),                                                                   (3*m2 - 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2))/(3*m4)];
[             m2,                    0, (m2*(3*m2 + 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2)))/(3*A1*m1*m4) - (7*m2^2 - 8*A1*m1*m3 + 2*A1*m1*m4*qcon)/(6*A1*m1*m4), (m2*(3*m2 - 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2)))/(3*A1*m1*m4) - (7*m2^2 - 8*A1*m1*m3 + 2*A1*m1*m4*qcon)/(6*A1*m1*m4)];
[              0,                   m2,                                                                                                                                       1,                                                                                                                                       1]
];

EVL = [
[                                                                               (u*(7*u^2 - 8*A1*Es + 14*A1*Y*qcon))/(8*A1*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),                                                                          -(3*(u^2 + 2*A1*Y*qcon))/(4*A1*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon)), (3*(u^2 + 2*A1*Y*qcon))/(4*rho*u*(u^2 - 2*A1*Es + 2*A1*Y*qcon)), (qcon*(u^2 - 8*A1*Es + 2*A1*Y*qcon))/(4*rho*u*(u^2 - 2*A1*Es + 2*A1*Y*qcon))];
[                                                                                                               (3*Y*u)/(4*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),                                                                                               -(3*Y)/(2*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),                (3*A1*Y)/(2*rho*u*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),        (2*u^2 - 4*A1*Es + A1*Y*qcon)/(2*rho*u*(u^2 - 2*A1*Es + 2*A1*Y*qcon))];
[  (3*2^(1/2)*Y*u*(4*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon) + 2^(1/2)*rho*u*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(1/2)))/(16*rho*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(3/2)), -(3*2^(1/2)*Y*(2*u^2 - 4*A1*Es + 4*A1*Y*qcon + 2^(1/2)*u*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(1/2)))/(8*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(3/2)),                     -(3*A1*Y)/(4*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),                              (3*A1*Y*qcon)/(4*(u^2 - 2*A1*Es + 2*A1*Y*qcon))];
[ -(3*2^(1/2)*Y*u*(4*rho*(u^2 - 2*A1*Es + 2*A1*Y*qcon) - 2^(1/2)*rho*u*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(1/2)))/(16*rho*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(3/2)),  (3*2^(1/2)*Y*(2*u^2 - 4*A1*Es + 4*A1*Y*qcon - 2^(1/2)*u*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(1/2)))/(8*(- u^2 + 2*A1*Es - 2*A1*Y*qcon)^(3/2)),                     -(3*A1*Y)/(4*(u^2 - 2*A1*Es + 2*A1*Y*qcon)),                              (3*A1*Y*qcon)/(4*(u^2 - 2*A1*Es + 2*A1*Y*qcon))]
];

EVL*EVR-diag([1,1,1,1])