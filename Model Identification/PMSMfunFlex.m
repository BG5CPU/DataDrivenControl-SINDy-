function [dxdt] = PMSMfunFlex(t, x, u, parameter)

% parameter = [A1,A2,A3,A4,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,D6,D7,D8];
A1 = parameter(1);
A2 = parameter(2);
A3 = parameter(3);
A4 = parameter(4);
C1 = parameter(5);
C2 = parameter(6);
C3 = parameter(7);
C4 = parameter(8);
C5 = parameter(9);
D1 = parameter(10);
D2 = parameter(11);
D3 = parameter(12);
D4 = parameter(13);
D5 = parameter(14);
D6 = parameter(15);
D7 = parameter(16);
D8 = parameter(17);

dx1 = x(2);
dx2 = x(3);
dx3 = A1*x(1)+A2*x(2)+A3*x(3)+A4*x(5);
dx4 = C1*x(4)+C2*x(1)*x(5)+C3*x(2)*x(5)+C4*x(3)*x(5)+C5*u(1);
dx5 = D1*x(5)+D2*x(1)*x(4)+D3*x(2)*x(4)+D4*x(3)*x(4)+D5*x(1)+D6*x(2)+D7*x(3)+D8*u(2);

dxdt = [dx1;dx2;dx3;dx4;dx5];

end











