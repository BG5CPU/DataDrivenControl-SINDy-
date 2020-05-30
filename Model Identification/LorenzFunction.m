function [dx] = LorenzFunction(t, x, parameter, u)

GAMA = parameter(1); 
SIGMA = parameter(2);

dx = [
    -x(1)+x(2)*x(3)+u(1);
    -x(2)-x(1)*x(3)+GAMA*x(3)+u(2);
    SIGMA*(x(2)-x(3))+u(3);
];
end

