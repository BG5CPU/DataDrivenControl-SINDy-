clear; clc;

%% Simulate the original system
R = 0.011;
L = 0.018;
Fi = 0.055;
p = 4;
JM = 0.002; JL = 0.1;
BM = 0.18; BS = 0.00;
KS = 500;

A1 = -(BM*KS)/(JM*JL);
A2 = -(JM*KS+JL*KS+BM*BS)/(JM*JL);
A3 = -(JL*BS+JM*BS+JL*BM)/(JM*JL);
A4 = (p*Fi)/(JM*JL);
C1 = -R/L; C2 = KS*p; C3 = BS*p; C4 = JL*p; C5 = 1/L;
D1 = -R/L; D2 = -KS*p; D3 = -BS*p; D4 = -JL*p;
D5 = (Fi*KS*p)/L; D6 = (Fi*BS*p)/L;
D7 = (Fi*JL*p)/L; D8 = 1/L;
parameter = [A1,A2,A3,A4,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,D6,D7,D8];
B3 = JL; B2 = BS; B1 = KS;
%initial condition
ic = [0, 0, 0, 1, 1]; 
% u(1) = ud
% u(2) = uq
u = [0, 0];
dt = 0.001;
tLimit = 100;
tspan = dt:dt:tLimit;
% tud = tspan;
% ud = 1*sin(2*pi*tud/20);

% Intergration
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, x] = ode45(@(t, x)PMSMfunFlex(t, x, u, parameter), tspan, ic, opts);

omega = B1*x(:,1)+B2*x(:,2)+B3*x(:,3); % this variable can be measured
id = x(:,4); % this variable can be measured
iq = x(:,5); % this variable can be measured
% plot3(id, iq, omega);
% plot(omega);

% creat Hankel matrix
rowNum = 10;  % changing the number of shift-stacked rows
% r = 3;  % number of modes to keep
nVars = 3; % number of modes to keep
% X = zeros(stackmax,size(x,1)-stackmax);
colNum = length(omega)-rowNum;
Hankel = zeros(rowNum, colNum);
for k=1:rowNum
    Hankel(k,:) = omega(k+1:k+colNum);
end
[U,S,V] = svd(Hankel,'econ');

% compute derivative using fourth order central difference
% use TD if more error 
for i=3:length(V)-3
    for k=1:nVars
        dV(i-2,k) = (1/(12*dt))*(-V(i+2,k)+8*V(i+1,k)-8*V(i-1,k)+V(i-2,k));
    end
end  
% concatenate
xx = V(3:end-3,1:nVars);
% xx = [xx,x(6:end-10,5)];
dxx = dV;

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 2; usesine = 0; nVars = 3;
Theta = poolData(xx, nVars, polyorder, usesine);
Theta = [Theta,x(6:end-10,5)];


% normalize columns of Theta (required in new time-delay coords)
mTheta = size(Theta,2);
normTheta = zeros(mTheta,1);
for k=1:mTheta
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
% compute Sparse regression: sequential least squares
% requires different lambda parameters for each column
clear Xi;
Xi(:,1) = sparsifyDynamics(Theta,dxx(:,1),0.8,1);
Xi(:,2) = sparsifyDynamics(Theta,dxx(:,2),20,1);
Xi(:,3) = sparsifyDynamics(Theta,dxx(:,3),50,1); 
% Theta = poolData(x,r,polyorder,usesine);
Xi
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end


% dt = 0.001;
% tLimit = 10;
% tspan = dt:dt:tLimit;
u_t = tspan;
u_q = iq;
ic = [0, 0, 1]; 
[tB,xB]=ode45(@(t,x)LoadFunc1(t, x, u_t, u_q),tspan,ic,opts);  % approximate

figure; plot3(xB(1000:end,1),xB(1000:end,2),xB(1000:end,3));

figure; plot3(V(1000:end,1),V(1000:end,2),V(1000:end,3));


















