function [k_LQ]=LQ(Cz,Dzu,A_a,B_a,t)
format long
nx=size(A_a,1);
nu=size(B_a,2);
nz=size(Cz,1);
if t==0
Q=Cz'*Cz;
R=Dzu'*Dzu;
P = sdpvar(nx,nx);
Y = sdpvar(nu,nx);
gamma = sdpvar(5,5);
F1=([P*A_a'+A_a*P+Y'*B_a'+B_a*Y P Y';
        P -inv(Q) zeros(5,1);
        Y zeros(1,5) -inv(R)] <= 0);
F2=(P >= 0);
F3=([gamma eye(5);eye(5) P] >= 0);
F=F1+F2+F3;
opts=sdpsettings('solver','sedumi','verbose',0);
solvesdp(F,trace(gamma),opts);
k_LQ=double(Y)*inv(double(P));
end
end