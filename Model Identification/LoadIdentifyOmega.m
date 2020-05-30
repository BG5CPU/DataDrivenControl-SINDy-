clear; clc;

%% Simulate the original system
R = 0.011;
L = 0.018;
Fi = 0.055;
p = 4;
JM = 0.002; JL = 0.1;
BM = 0.18; BS = 0.002;
KS = 500;

a1 = -(BM*KS)/(JM*JL);
a2 = -(JM*KS+JL*KS+BM*BS)/(JM*JL);
a3 = -(JL*BS+JM*BS+JL*BM)/(JM*JL);
a4 = (p*Fi)/(JM*JL);
c3 = JL; c2 = BS; c1 = KS;

A1 = -c1/c3; A2 = -c2/c3; A3 = 1/c3;
B1 = a1*c3-a3*c1-c1*c2/c3;
B2 = a2*c3-a3*c2+c1-c2*c2/c3;
B3 = a3+c2/c3; B4 = a4*c3;
C1 = -R/L; C2 = 1*p; C3 = 1/L;
D1 = -R/L; D2 = -1*p; D3 = Fi*p/L; D4 = 1/L;
parameter = [A1; A2; A3; B1; B2; B3; B4; C1; C2; C3; D1; D2; D3; D4];
%initial condition
ic = [0, 0, 1, 1, 1]; 
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
[t, x] = ode45(@(t, x)PMSMfunFlexXX(t, x, u, parameter), tspan, ic, opts);

% figure; plot3(x(1000:end,1),x(1000:end,2),x(1000:end,3));

omega = x(:,3); % this variable can be measured
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
% figure; plot3(V(1000:end,1),V(1000:end,2),V(1000:end,3));

VV = [];
for i=1:length(V)
    if mod(i,10) == 1
        VV = [VV; V(i,:)];
    end
end
figure('color',[1 1 1]);
set(gcf,'position',[50 50 800 300]);
title('Original System','Fontname', 'Times New Roman','FontSize',15);
subplot(1,3,1);
plot(VV(1000:end,1),VV(1000:end,2),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\bar{x}_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\bar{x}_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.06 0.18 0.25 0.75]);
subplot(1,3,2);
plot(VV(1000:end,1),VV(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\bar{x}_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\bar{x}_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.39 0.18 0.25 0.75]);
subplot(1,3,3);
plot(VV(1000:end,2),VV(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\bar{x}_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\bar{x}_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.72 0.18 0.25 0.75]);

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
polyorder = 3; usesine = 0; nVars = 3;
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
Xi(:,1) = sparsifyDynamics(Theta,dxx(:,1),1,1);
Xi(:,2) = sparsifyDynamics(Theta,dxx(:,2),40,1);
Xi(:,3) = sparsifyDynamics(Theta,dxx(:,3),10,1); 
% Theta = poolData(x,r,polyorder,usesine);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end
Xi

% dt = 0.001;
% tLimit = 10;
% tspan = dt:dt:tLimit;
u_t = tspan;
u_q = iq;
ic = [0, 0, 1]; 
[tB,xB]=ode45(@(t,x)LoadFunc1(t, x, u_t, u_q),tspan,ic,opts);  % approximate

% figure; plot3(xB(1000:end,1),xB(1000:end,2),xB(1000:end,3));
xxB = [];
for i=1:length(V)
    if mod(i,10) == 1
        xxB = [xxB; xB(i,:)];
    end
end
figure('color',[1 1 1]);
set(gcf,'position',[50 50 800 300]);
title('Original System','Fontname', 'Times New Roman','FontSize',15);
subplot(1,3,1);
plot(xxB(1000:end,1),xxB(1000:end,2),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\tilde{x}_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\tilde{x}_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.06 0.18 0.25 0.75]);
subplot(1,3,2);
plot(xxB(1000:end,1),xxB(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\tilde{x}_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\tilde{x}_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.39 0.18 0.25 0.75]);
subplot(1,3,3);
plot(xxB(1000:end,2),xxB(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$\tilde{x}_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$\tilde{x}_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.72 0.18 0.25 0.75]);



















