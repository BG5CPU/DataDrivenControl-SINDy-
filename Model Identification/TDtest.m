clear; clc;

%% generate Data
polyorder = 3;
usesine = 0;
nVars = 3;

% system parameters
Rs = 0.011;
Ld = 0.018;
Lq = 0.018;
Fi = 0.055;
p = 4;
J = 0.102;
b = 0.18;
tao = Lq/Rs;
k = b/(tao*p*p*Fi);
GAMA =  Fi/(k*Lq)
SIGMA = tao*b/J
parameter = [GAMA SIGMA];

ic = [1, 1, 1]; %initial condition
% u(1) = ud
% u(2) = uq
% u(3) = TL
u = [0; 0; 0];
dt = 0.001;
tLimit = 100;
tspan = dt:dt:tLimit;

% Integrate
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, x] = ode45(@(t, x)LorenzFunction(t, x, parameter, u), tspan, ic, opts);
xclean = x;
% add noise
eps = 10;
x = x + eps*randn(size(x));

% TD filter
r = 100000; h = dt; h0 = 50*h; h1 = 1.2*h0;
[Xtrack, Xdiff] = TDfunc(x, r, h, h0, h1);
% XXtrack = Xtrack(10000:90000,:);
% XXdiff = Xdiff(10000:90000,:);

% compute clean derivative  (just for comparison!)
dxclean = zeros(size(xclean));
for i=1:length(xclean)
    dxclean(i,:) = LorenzFunction(0, xclean(i,:), parameter, u);
end

% 
xt(:,1) = cumsum(Xdiff(:,1))*dt;
xt(:,2) = cumsum(Xdiff(:,2))*dt;
xt(:,3) = cumsum(Xdiff(:,3))*dt;
xt(:,1) = xt(:,1) - (mean(xt(2000:end-2000,1)) - mean(x(2000:end-2000,1)));
xt(:,2) = xt(:,2) - (mean(xt(2000:end-2000,2)) - mean(x(2000:end-2000,2)));
xt(:,3) = xt(:,3) - (mean(xt(2000:end-2000,3)) - mean(x(2000:end-2000,3)));
xt = xt(2000:end-2000,:);
dxt = Xdiff(2000:end-2000,:);  % trim off ends (overly conservative)

% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(xt, nVars, polyorder, usesine);
m = size(Theta,2);

% compute Sparse regression: sequential least squares
lambda = 0.3;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta, dxt, lambda, nVars)

% 
tspan = 0.005:0.005:10;
[tA,xA] = ode45(@(t, x)LorenzFunction(t, x, parameter, u), tspan, ic, opts);   % true model
[tB,xB] = ode45(@(t,x)sparseGalerkin(t, x, Xi, polyorder, usesine), tspan, ic, opts);  % approximate

xDelt = abs(xB(:,3)-xA(:,3));

%% figure
% figure;
% plot(t, x(:,1), '-', 'LineWidth', 0.5, 'color', [0, 1, 0]);
% hold on;
% plot(t, Xtrack(:,1), '-', 'LineWidth', 0.5, 'color', [1, 0, 0]);
% hold on;
% plot(t, xclean(:,1), '-', 'LineWidth', 0.5, 'color', [0, 0, 1]);
% 
% figure;
% plot(t, Xdiff(:,1), '-', 'LineWidth', 0.5, 'color', [0, 0, 1]);
% hold on;
% plot(t, dxclean(:,1), '-', 'LineWidth', 0.5, 'color', [1, 0, 0]);

xx1 = xB(floor(end/2):end,1);
xx2 = xB(floor(end/2):end,2);
xx3 = xB(floor(end/2):end,3);

% figure('color',[1 1 1]);
% set(gcf,'position',[50 50 500 400]);
% plot3(xx1, xx2, xx3, '-', 'LineWidth', 0.75, 'color', [0, 0, 0]);
% set(gca,'fontsize',15,'fontname','Times');
% title('Identified System, $\eta=0.1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% % title('Original System','interpreter','latex','Fontname', 'Times New Roman','FontSize',16)
% xlabel('$x_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16,'FontAngle','italic');
% ylabel('$x_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16,'FontAngle','italic');
% zlabel('$x_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16,'FontAngle','italic');
% grid on;
% set(gca,'position',[0.14 0.10 0.80 0.82]);

figure('color',[1 1 1]);
set(gcf,'position',[50 50 500 600]);
subplot(3,1,1); % first
plot(tA, xA(:,1), ':', 'LineWidth', 2.5, 'color', [1, 0, 0]);
hold on;
plot(tB, xB(:,1), '-', 'LineWidth', 1, 'color', [0, 0, 0]);
set(gca,'fontsize',15,'fontname','Times');
ylabel('$x_1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16);
grid on;
set(gca,'position',[0.11 0.70 0.87 0.25]);
title('$\eta=0.1$','interpreter','latex','Fontname', 'Times New Roman','FontSize',15);
subplot(3,1,2); % second
plot(tA, xA(:,2), ':', 'LineWidth', 2.5, 'color', [1, 0, 0]);
hold on;
plot(tB, xB(:,2), '-', 'LineWidth', 1, 'color', [0, 0, 0]);
set(gca,'fontsize',15,'fontname','Times');
ylabel('$x_2$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16);
grid on;
set(gca,'position',[0.11 0.40 0.87 0.25]);
legend('Original System', 'Identified System', 'Fontname', 'Times New Roman','FontSize',16);
subplot(3,1,3); % third
plot(tA, xA(:,3), ':', 'LineWidth', 2.5, 'color', [1, 0, 0]);
hold on;
plot(tB, xB(:,3), '-', 'LineWidth', 1, 'color', [0, 0, 0]);
set(gca,'fontsize',15,'fontname','Times');
ylabel('$x_3$','interpreter','latex','Fontname', 'Times New Roman','FontSize',16);
grid on;
set(gca,'position',[0.11 0.10 0.87 0.25]);
xlabel('time(s)','interpreter','latex','Fontname', 'Times New Roman','FontSize',16);


figure('color',[1 1 1]);
set(gcf,'position',[50 50 800 300]);
% xDelt = abs(xB(:,3)-xA(:,3));
load eta10.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [0.098, 0.098, 0.439]);
hold on;
load eta8.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [0, 0, 0.8039]);
hold on;
load eta6.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [0.5, 0.3750, 1]);
hold on;
load eta4.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [0.5, 0.5, 0.5]);
hold on;
load eta2.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [1, 0, 0]);
hold on;
load eta1.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [1, 1, 0]);
hold on;
load eta0.1.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [1, 0.64706, 0]);
hold on;
load eta0.01.mat
semilogy(tspan, xDelt, '-', 'LineWidth', 1, 'color', [1, 0.411, 0.705]);
hold on;
legend('$\eta=10$', '$\eta=8$', '$\eta=6$', '$\eta=4$', '$\eta=2$', '$\eta=1$', '$\eta=0.1$', '$\eta=0.01$', 'interpreter','latex','Fontname', 'Times New Roman','FontSize',12);
grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('time(s)','interpreter','latex','Fontname', 'Times New Roman','FontSize',12);
ylabel('Error','interpreter','latex','Fontname', 'Times New Roman','FontSize',12);
set(gca,'position',[0.078 0.14 0.78 0.82]);


