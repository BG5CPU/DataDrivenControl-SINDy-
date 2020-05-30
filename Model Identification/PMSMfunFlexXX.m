function [dxdt] = PMSMfunFlexXX(t, x, u, parameter)

% parameter = [A1; A2; A3; B1; B2; B3; B4; C1; C2; C3; D1; D2; D3; D4];
A1 = parameter(1);
A2 = parameter(2);
A3 = parameter(3);
B1 = parameter(4);
B2 = parameter(5);
B3 = parameter(6);
B4 = parameter(7);
C1 = parameter(8);
C2 = parameter(9);
C3 = parameter(10);
D1 = parameter(11);
D2 = parameter(12);
D3 = parameter(13);
D4 = parameter(14);


dx1 = x(2);
dx2 = A1*x(1)+A2*x(2)+A3*x(3);
dx3 = B1*x(1)+B2*x(2)+B3*x(3)+B4*x(5);
dx4 = C1*x(4)+C2*x(3)*x(5)+C3*u(1);
dx5 = D1*x(5)+D2*x(3)*x(4)+D3*x(3)+D4*u(2);

dxdt = [dx1;dx2;dx3;dx4;dx5];

end











