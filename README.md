Data-Driven Control for Permanent Magnet Synchronous Motor(PMSM) and PMSM with flexible load

In "Model Identification", Sparse Identification of Nonlinear Dynamical Systems(SINDy) is used to identify the model parameters of PMSM. The identified parameters can be used to control.

In "State Feedback Control", 

...\State Feedback Control\controlMotorSelfK3K4K5K6.m  is used to calculate the state feedback gain K and error feedback matrix G(in state observer) according to the desired position of the poles.

...\State Feedback Control\*.slx is the state feedback control for the PMSM and PMSM with flexible load.


About the SINDy, 
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% https://doi.org/10.1073/pnas.1517384113
% http://www.databookuw.com/
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz