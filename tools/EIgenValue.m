% matlab ���룺�����������
% ����
syms rho u Es Y real
syms m_1 m_2 m_3 m_4 real
syms A1 qcon real

% ��ֵ
P = m_3-m_2*m_2/m_1/2-qcon*m_4;
Um = [m_1 m_2 m_3 m_4];
Fm = [m_2; m_2*m_2/m_1+P/3; m_3*m_2/m_1+m_2*P/m_1/3; m_2*m_4/m_1];

% ���
J = jacobian(Fm,Um);
[RM,RD] = eig(J);
[LM,LD] = eig(J');
 
% EVR = [
% [ (2*A1*m1^2)/m2, -(2*A1*m1^2*qcon)/m2,                                                                                                                                   m1/m4,                                                                                                                                   m1/m4];
% [     (2*A1*m1) ,     -(2*A1*m1*qcon) ,                                                                   (3*m2 + 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2))/(3*m4),                                                                   (3*m2 - 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2))/(3*m4)];
% [             m2,                    0, (m2*(3*m2 + 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2)))/(3*A1*m1*m4) - (7*m2^2 - 8*A1*m1*m3 + 2*A1*m1*m4*qcon)/(6*A1*m1*m4), (m2*(3*m2 - 2^(1/2)*(- m2^2 + 2*A1*m1*m3 - 2*A1*m1*m4*qcon)^(1/2)))/(3*A1*m1*m4) - (7*m2^2 - 8*A1*m1*m3 + 2*A1*m1*m4*qcon)/(6*A1*m1*m4)];
% [              0,                   m2,                                                                                                                                       1,                                                                                                                                       1]
% ];

% % ��֤һ��|A-��I|����ʽֵΪ0 
% for i=1:4
%     det(J-RD(i,i)*diag([1,1,1,1],0))
% end
 
% ��֤����Ax-��x=0 && yA-��x=0
for i=1:4
    rr = J*RM(:,i) -RD(i,i).*RM(:,i) ;simplify(rr)
%     ll = LM(:,i)'*J-RD(i,i).*LM(:,i)';simplify(ll)
end

%% ֤�� RM��������������Ϊ����������LM��������������ת��Ϊ��������