function Km = calc_Kmin(v)
%%

m = 1575;
Iz = 2875;
Lf = 1.2;
Lr = 1.6;
Cf = 19000;
Cr = 33000;

% Specify a state-space model, |G(s)|, of the lateral vehicle dynamics.
A = [   0               1             0           0
        0    -(2*Cf+2*Cr)/m/v        0    -v-(2*Cf*Lf-2*Cr*Lr)/m/v;...
        0               0             0           1; ...
        0   -(2*Cf*Lf-2*Cr*Lr)/Iz/v  0    -(2*Cf*Lf^2+2*Cr*Lr^2)/Iz/v];
B = [0  2*Cf/m 0 2*Cf*Lf/Iz]';
C = eye(4);
D = zeros(4,1);
%CreaCurvatura;
load('BusSignals2.mat')
% load('qw.mat');


%% H8 statico
Bu = B;
Bw = B;
alfa1 = 10^0;
Cz = [0  0 0 0  alfa1*1;
      0  0 0 0 0;];

rho=10^-1;
Dzu = rho * [0; 1];
Dzw =  zeros(size(Cz,1),2);

% sistema aumentato
[n,m] = size(A);
new_row1 = [0 0 -C(3,3) 0 0];

A_a = [A zeros(n,1);new_row1];
[n,m] = size(B);
Bu_a = [Bu; zeros(1,m)];
[n,m] = size(C);
C_a = [C zeros(n,1);
       zeros(1,m) eye(1)];
[n,m] = size(D);
D_a = [D; zeros(1,m)];
[n,m] = size(Bw);
Bw_a = [Bw zeros(n,1);
        zeros(1,m) eye(1,1)];

%
Km = h8(A_a, Bu_a, Bw_a, Cz, Dzu, Dzw)
end
