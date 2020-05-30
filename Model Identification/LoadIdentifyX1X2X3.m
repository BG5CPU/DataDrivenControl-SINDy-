clear; clc;

%% Simulate the original system
R = 0.011;
L = 0.018;
Fi = 0.055;
p = 4;
JM = 0.002; JL = 0.1;
BM = 0.18; BS = 0.002;
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
dt = 0.01;
tLimit = 100;
tspan = dt:dt:tLimit;
% tud = tspan;
% ud = 1*sin(2*pi*tud/20);

% Intergration
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, x] = ode45(@(t, x)PMSMfunFlex(t, x, u, parameter), tspan, ic, opts);

id = x(:,4); % this variable can be measured
iq = x(:,5); % this variable can be measured
omega = B1*x(:,1)+B2*x(:,2)+B3*x(:,3); % this variable can be measured
% plot(omega);
% figure; plot3(x(1000:end,1),x(1000:end,2),x(1000:end,3));
figure('color',[1 1 1]);
set(gcf,'position',[50 50 800 300]);
title('Original System','Fontname', 'Times New Roman','FontSize',15);
subplot(1,3,1);
plot(x(1000:end,1),x(1000:end,2),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$x_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$x_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.06 0.18 0.25 0.75]);
subplot(1,3,2);
plot(x(1000:end,1),x(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$x_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$x_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.39 0.18 0.25 0.75]);
subplot(1,3,3);
plot(x(1000:end,2),x(1000:end,3),'-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
set(gca,'fontsize',12,'fontname','Times');
xlabel('$x_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
ylabel('$x_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
grid on; set(gca,'position',[0.72 0.18 0.25 0.75]);

xclean = x;
dxclean = zeros(size(xclean));
for i=1:length(xclean)
    dxclean(i,:) = PMSMfunFlex(t, xclean(i,:), u, parameter);
end
xclean = xclean(:,1:3);
dxclean = dxclean(:,1:3);

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 3; usesine = 0; nVars = 3;
Theta = poolData(xclean, nVars, polyorder, usesine);
Theta = [Theta, x(:,5)];

% normalize columns of Theta (required in new time-delay coords)
mTheta = size(Theta,2);
for k=1:mTheta
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
% compute Sparse regression: sequential least squares
% requires different lambda parameters for each column
Xi(:,1) = sparsifyDynamics(Theta,dxclean(:,1),0.5,1);
Xi(:,2) = sparsifyDynamics(Theta,dxclean(:,2),0.5,1);
Xi(:,3) = sparsifyDynamics(Theta,dxclean(:,3),50,1); 
% Theta = poolData(x,r,polyorder,usesine);
for k=1:length(Xi)
    Xi(k,:) = Xi(k,:)/normTheta(k);
end
Xi






















