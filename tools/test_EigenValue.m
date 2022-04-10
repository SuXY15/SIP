%% Get the eigenvalues and eigenvectors by symbolic computation
clc;
% definition
syms m1 m2 m3 m4 real
syms Q gamma real
syms m1 m3 m4 Q gamma positive

% declaration
rho = m1;
u = m2/m1;
E = m3/m1;
Y = m4/m1;
e = (E - 1/2*u*u - Q*Y);
p = (gamma-1)*rho*e;
Um = [m1, m2, m3, m4]; % [rho; rho*u; rho*E; rho*Y]
Fm = [rho*u; rho*u*u+p; rho*u*E+u*p; rho*u*Y];

% solving
J = jacobian(Fm,Um);
[RM,RD] = eig(J);
[LM,LD] = eig(J');

% validation for eigenvalues£º|A-ŠËI|=0
for i=1:4
    v = det(J-RD(i,i)*diag([1,1,1,1],0)); disp(v');
end

% validation for eigenvectors£ºAx-ŠËx=0 && yA-ŠËx=0
for i=1:4
    rr = J*RM(:,i)-RD(i,i).*RM(:,i); v = simplify(rr); disp(v');
    ll = LM(:,i)'*J-LD(i,i).*LM(:,i)'; v = simplify(ll); disp(v);
end

%% get simplified formula by symbolic operation
clc; clear;
syms u real
syms rho E Y real positive
syms Q gamma real positive
m1 = rho;
m2 = rho*u;
m3 = rho*E;
m4 = rho*Y;
e = (E - 1/2*u*u - Q*Y);
p = (gamma-1)*rho*e;
c = sqrt(gamma*p/rho);
H = E + p/rho;
gamma = 4/3;

J = simplify([
  [                                                                                                                0,                                                                                 1,                           0,                      0];
  [ m1*(gamma - 1)*(m2^2/m1^3 - m3/m1^2 + (Q*m4)/m1^2) - m2^2/m1^2 - (gamma - 1)*(m2^2/(2*m1^2) - m3/m1 + (Q*m4)/m1),                                                   (2*m2)/m1 - (m2*(gamma - 1))/m1,                   gamma - 1,         -Q*(gamma - 1)];
  [                                                m2*(gamma - 1)*(m2^2/m1^3 - m3/m1^2 + (Q*m4)/m1^2) - (m2*m3)/m1^2, m3/m1 - (gamma - 1)*(m2^2/(2*m1^2) - m3/m1 + (Q*m4)/m1) - (m2^2*(gamma - 1))/m1^2, m2/m1 + (m2*(gamma - 1))/m1, -(Q*m2*(gamma - 1))/m1];
  [                                                                                                    -(m2*m4)/m1^2,                                                                             m4/m1,                           0,                  m2/m1];
])

RD = simplify([
  [ m2/m1,     0,                                                                                     0,                                                                                     0];
  [     0, m2/m1,                                                                                     0,                                                                                     0];
  [     0,     0, (2*m2 + 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*m1),                                                                                     0];
  [     0,     0,                                                                                     0, (2*m2 - 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*m1)];
])

RM = simplify([
  [ (2*m1^2)/m2^2, -(2*Q*m1^2)/m2^2,                                                                                                                                                                       m1/m4,                                                                                                                                                                       m1/m4];
  [     (2*m1)/m2,     -(2*Q*m1)/m2,                                                                                       (2*m2 + 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*m4),                                                                                       (2*m2 - 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*m4)];
  [             1,                0, (m2*(2*m2 + 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2)))/(2*m1*m4) - (gamma*m2^2 + m2^2 - 2*Q*m1*m4 - 2*gamma*m1*m3 + 2*Q*gamma*m1*m4)/(2*m1*m4), (m2*(2*m2 - 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2)))/(2*m1*m4) - (gamma*m2^2 + m2^2 - 2*Q*m1*m4 - 2*gamma*m1*m3 + 2*Q*gamma*m1*m4)/(2*m1*m4)];
  [             0,                1,                                                                                                                                                                           1,                                                                                                                                                                           1];
])

LM = simplify([
  [ (gamma*m2^2 + m2^2 - 2*Q*m1*m4 - 2*gamma*m1*m3 + 2*Q*gamma*m1*m4)/(2*m1^2), -m4/m1, (m2*(2*m2 + 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2)))/(2*Q*m1^2*(gamma - 1)) - (gamma*m2^2 + m2^2)/(2*Q*m1^2*(gamma - 1)), (m2*(2*m2 - 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2)))/(2*Q*m1^2*(gamma - 1)) - (gamma*m2^2 + m2^2)/(2*Q*m1^2*(gamma - 1))];
  [                                                                     -m2/m1,      0,                     (gamma*m2)/(Q*m1*(gamma - 1)) - (2*m2 + 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*Q*m1*(gamma - 1)),                     (gamma*m2)/(Q*m1*(gamma - 1)) - (2*m2 - 2^(1/2)*gamma^(1/2)*(-(gamma - 1)*(m2^2 - 2*m1*m3 + 2*Q*m1*m4))^(1/2))/(2*Q*m1*(gamma - 1))];
  [                                                                          1,      0,                                                                                                                                                    -1/Q,                                                                                                                                                    -1/Q];
  [                                                                          0,      1,                                                                                                                                                       1,                                                                                                                                                       1];
])

%% get the results for eigenvalues and eigenvectors
clc; clear;
syms u real
syms rho E Y real positive
syms Q gamma real positive
% rho = 2;
% u = 3;
% E = 100;
% Y = 1;
% gamma = 4/3;
% Q = 0.1;

e = (E - 1/2*u*u - Q*Y);
p = (gamma-1)*rho*e;
c = sqrt(gamma*p/rho);
H = E + p/rho;

K1 = (gamma-1)*u*u/2;
K2 = (gamma-1)*(-u);
K3 = (gamma-1);
K4 = (gamma-1)*(-Q);

J = [
  [        0,      1,         0,    0];
  [  -u*u+K1, 2*u+K2,        K3,   K4]
  [-u*H+u*K1, H+u*K2,  u*(1+K3), u*K4]
  [     -u*Y,      Y,         0,    u];
];

D = [
  [u-c, 0, 0,   0];
  [  0, u, 0,   0];
  [  0, 0, u,   0];
  [  0, 0, 0, u+c];
];

% [V_, D_, W_] = eig(J);

RM = [
  [1,     1,      0,     1];
  [u-c,   u,      0,   u+c];
  [H-u*c, u*u/2,  Q, H+u*c];
  [Y,     0,      1,     Y];
];

c2 = 2*c;
c22 = 2*c*c;
LM = [
  [u/c2+K1/c22,   -1/c2-K3*u/c22,   K3/c22,    K4/c22];
  [1-K1/c/c,       K3*u/c/c,       -K3/c/c,   -K4/c/c];
  [-Y*K1/c/c,    Y*K3*u/c/c,     -Y*K3/c/c, 1-Y*K4/c/c];
  [-u/c2+K1/c22,   1/c2-K3*u/c22,   K3/c22,    K4/c22];
]';

% % validation for eigenvalues£º|A-ŠËI|=0
% for i=1:4
%     v = det(J-RD(i,i)*diag([1,1,1,1],0)); disp(v');
% end

% validation for eigenvectors£ºAx-ŠËx=0 && yA-ŠËx=0
for i=1:4
    rr = J*RM(:,i)-D(i,i).*RM(:,i); v = simplify(rr); disp(v');
    ll = LM(:,i)'*J-D(i,i).*LM(:,i)'; v = simplify(ll); disp(v);
end
simplify(LM'*RM)

%% get the final eigen-space RF and RF^-1
clc; clear;
syms u real
syms rho E Y real positive
syms Q gamma real positive

e = (E - 1/2*u*u - Q*Y);
p = (gamma-1)*rho*e;
c = sqrt(gamma*p/rho);
H = E + p/rho;

K1 = (gamma-1)*u*u/2;
K2 = (gamma-1)*(-u);
K3 = (gamma-1);
K4 = (gamma-1)*(-Q);

J = [
  [        0,      1,         0,    0];
  [  -u*u+K1, 2*u+K2,        K3,   K4]
  [-u*H+u*K1, H+u*K2,  u*(1+K3), u*K4]
  [     -u*Y,      Y,         0,    u];
];

D = [
  [u-c, 0, 0,   0];
  [  0, u, 0,   0];
  [  0, 0, u,   0];
  [  0, 0, 0, u+c];
];

RF = [
  [1/c,       1/c,      0,     1/c];
  [u/c-1,     u/c,      0,   u/c+1];
  [H/c-u, u*u/2/c,    Q/c,   H/c+u];
  [Y/c,         0,    1/c,     Y/c];
];

g1 = gamma-1;
c2 = 2*c;
c4 = 4*c;
LF = [
  [ u/2+g1*u*u/c4,   -1/2-g1*u/c2,   g1/c2,    -g1*Q/c2];
  [   c-g1*u*u/c2,        g1*u/c,   -g1/c,      g1*Q/c ];
  [  -Y*g1*u*u/c2,      Y*g1*u/c, -Y*g1/c,  c+Y*g1*Q/c ];
  [-u/2+g1*u*u/c4,    1/2-g1*u/c2,   g1/c2,    -g1*Q/c2];
];

simplify(RF*LF)        % output is [ 1, 0, 0, 0]
                       %           [ 0, 1, 0, 0]
                       %           [ 0, 0, 1, 0]
                       %           [ 0, 0, 0, 1]

simplify(RF*D*LF - J)  % output is [ 0, 0, 0, 0]
                       %           [ 0, 0, 0, 0]
                       %           [ 0, 0, 0, 0]
                       %           [ 0, 0, 0, 0]
