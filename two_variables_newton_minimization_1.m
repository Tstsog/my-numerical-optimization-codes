% This Matlab code performs an optimization (finding the maximum/minimum) for two-variable function f(x1,x2) 
% using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The two-variable Rosenbrock's function: f(x1,x2) = 3*(1. - x1)^2*exp(-x1^2 - (x2+1)^2) - ...
%                                                   10*(0.5*x1 - x1^3 - x2^4)*exp(-x1^2 - x2^2) - ...
%                                                (1/3)*exp(-(x1+1)^2 - x2^2).   
%
% The Newton-Raphson iterative scheme: x^(k+1) = x^(k) - [Hessian(x^(k)]^(-1)*Jacobian(x^(k)), 
% where Jacobian, J(x) = df/dx, and Hessian = d^2f/dx^2 = d(Jacobian)/dx; 
% and 'k' defines the k-th iteration.
%
% The first and second order derivatives are taken with finite difference scheme.  
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% March 26, 2025 & University of North Dakota 
%
function [] = two_variables_newton_minimization_1
clc; clear two_variables_newton_minimization_1
%
format long 
%
x1 = -0.20; tol = 1e-6;
x2 = -1.00;
dx = 0.001; iter_max = 50.;
%

for iter = 1:iter_max
    %
    x_val = [x1;
             x2];
    %
    [iter, x1, x2, function_f(x1,x2)]
    %
    Jacobian_11 = (function_f(x1+dx,x2) - function_f(x1-dx,x2))/(2*dx); 
    Jacobian_22 = (function_f(x1,x2+dx) - function_f(x1,x2-dx))/(2*dx); 
    %
    Hessian_11 = (function_f(x1+dx,x2) - 2*function_f(x1,x2) + function_f(x1-dx,x2))/dx^2; 
    Hessian_12 = (function_f(x1+dx,x2+dx) - function_f(x1+dx,x2-dx) - function_f(x1-dx,x2+dx) + function_f(x1-dx,x2-dx))/(4*dx*dx);
    Hessian_21 = (function_f(x1+dx,x2+dx) - function_f(x1+dx,x2-dx) - function_f(x1-dx,x2+dx) + function_f(x1-dx,x2-dx))/(4*dx*dx); 
    Hessian_22 = (function_f(x1,x2+dx) - 2*function_f(x1,x2) + function_f(x1,x2-dx))/dx^2;     
    %
    Jabobian = [Jacobian_11;
                Jacobian_22];
    %
    Hessian = [Hessian_11, Hessian_12;
               Hessian_21, Hessian_22 ];    
    %
    x_val = x_val - Hessian\Jabobian; % x_n+1 = x_n - Jacobian/Hessian
    x1 = x_val(1);
    x2 = x_val(2);
    %
    if ((abs(Jacobian_11)) <= tol)
        break;
    end
%
end
%%%
%[iter, x1, x2, function_f(x1,x2)]
[
 1.000000000000000  -0.200000000000000  -1.000000000000000   7.945675490030952
 2.000000000000000  -0.350326843974732  -1.389875948197332   9.076280974322486
 3.000000000000000  -0.363782488079547  -1.259125943129218   9.276128416384399
 4.000000000000000  -0.365178904366060  -1.263317291714486   9.276396941728319
 5.000000000000000  -0.365185090079313  -1.263316126196002   9.276396942156921    
];

%%%
xx1 = -3.:0.1:3; xx2 = -3:0.1:3.;
[xx, yy] = meshgrid(xx1,xx2);

f_xy = 3.*(1. - xx).^2.*exp(-xx.^2 - (yy+1).^2) - ...
       10.*(0.5.*xx - xx.^3 - yy.^4).*exp(-xx.^2 - yy.^2) - ...
       (1./3).*exp(-(xx+1).^2 - yy.^2);

figure(1)
mesh(xx, yy, f_xy)
hold on
plot3(x1, x2, function_f(x1,x2),'r.','MarkerSize',15)
hold off
xlabel('$x_{1}$','interpreter','latex')
ylabel('$x_{2}$','interpreter','latex')
zlabel('$f(x_{1},x_{2})$','interpreter','latex', 'Rotation', 1)
set(gca,'FontSize',14)

%%%
return
end

%%%
%
function f = function_f(x1,x2)
%
f = 3*(1. - x1)^2*exp(-x1^2 - (x2+1)^2) - ...
    10*(0.5*x1 - x1^3 - x2^4)*exp(-x1^2 - x2^2) - ...
    (1/3)*exp(-(x1+1)^2 - x2^2);
%
return
end

