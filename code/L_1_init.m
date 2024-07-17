clear
clc
% To describe the lateral vehicle dynamics, this example uses a _bicycle
% model_ with the following parameters:
%
% * |m| is the total vehicle mass (kg).
% * |Iz| is the yaw moment of inertia of the vehicle (mNs^2).
% * |Lf| is the longitudinal distance from the center of gravity to the
% front tires (m).
% * |Lr| is the longitudinal distance from center of gravity to the rear
% tires (m).
% * |Cf| is the cornering stiffness of the front tires (N/rad).
% * |Cr| is the cornering stiffness of the rear tires (N/rad).
%
m = 1575;
Iz = 2875;
Lf = 1.2;
Lr = 1.6;
Cf = 19000;
Cr = 33000;
% Specify the longitudinal velocity in m/s.
Vx = 20;

%%
% Specify a state-space model, |G(s)|, of the lateral vehicle dynamics.
A = [   0               1             0           0
        0    -(2*Cf+2*Cr)/m/Vx        0    -Vx-(2*Cf*Lf-2*Cr*Lr)/m/Vx;...
        0               0             0           1; ...
        0   -(2*Cf*Lf-2*Cr*Lr)/Iz/Vx  0    -(2*Cf*Lf^2+2*Cr*Lr^2)/Iz/Vx];
B = [0  2*Cf/m 0 2*Cf*Lf/Iz]';
C = eye(4);
D = zeros(4,1);
x_0 = [0 0 0 0]';
%CreaCurvatura;
load('BusSignals2.mat')
% load('qw.mat');


%% L1 statico
Bu = B;
Bw = B;
alfa1 = 10^0;
Cz = [0  0 0 0  alfa1*1;
      0  0 0 0 0;];

rho=10^-10;
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
lambda=2.3972;

% k_1 gain
K=L_1(A_a,Bu_a,Bw_a,Cz,Dzu,Dzw,lambda)

Kfb = K(1,1:end-1)
Kff = K(1,5:end)