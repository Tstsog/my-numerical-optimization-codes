% This Matlab code performs an optimization (finding the maximum/minimum) for two-variable function f(x1,x2) 
% using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The two-variable Rosenbrock's function: f(x1,x2) = (1-x1)^2 + 100*(x2 - x1^2)^2.   
% The Newton-Raphson iterative scheme: x^(k+1) = x^(k) - [Hessian(x^(k)]^(-1)*Jacobian(x^(k)), 
% where Jacobian, J(x) = df/dx, and Hessian = d^2f/dx^2 = d(Jacobian)/dx; 
% and 'k' defines the k-th iteration.
%
% The first and second order derivatives are taken with finite difference scheme.  
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% March 21, 2025 & University of North Dakota 
%
function [] = two_variables_newton_minimization_Rosenbrock_function
clc; clear two_variables_newton_minimization_Rosenbrock_function; 
%
format long 
%
x1 = -1.50; tol = 1e-6;
x2 = 2.50;
dx = 0.001; iter_max = 50.;
%

for iter = 1:iter_max
    %
    x_val = [x1;
             x2];
    %
    output = [iter, x1, x2, function_f(x1,x2)]
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

% output = [iter, x1, x2, function_f(x1,x2)]
%1.000000000000000  -1.500000000000000   2.500000000000000  12.500000000000000
%2.000000000000000  -1.551026634742158   2.403079904226264   6.508414825621576
%3.000000000000000   0.1265535088818    -2.7982593476589   792.7773641608677
%4.000000000000000   0.128102525971115   0.016407860012014   0.760205205792582
%5.000000000000000   0.999469483721995   0.239658873831398  57.650669076590077
%6.000000000000000   0.999471646688536   0.998943572531390   0.000000279157222
%7.000000000000000   0.999800072827863   0.999600077762871   0.000000039972038
%8.000000000000000   0.999800039989454   0.999600119962913   0.000000039984006

%%%
xx1 = -1.5:0.1:2; xx2 = -1.5:0.1:2.;
[x1_plot, x2_plot] = meshgrid(xx1,xx2);
%
Rosenbloch_func = (1-x1_plot).^2 + 100.*(x2_plot - x1_plot.^2).^2;

figure(1)
mesh(x1_plot, x2_plot, Rosenbloch_func)
hold on
plot3(x1, x2, function_f(x1,x2),'r.','MarkerSize',15)
hold off
xlabel('$x_{1}$','interpreter','latex')
ylabel('$x_{2}$','interpreter','latex')
zlabel('$f(x_{1},x_{2})$','interpreter','latex', 'Rotation', 1)
set(gca,'FontSize',14)
%box on

%%%
return
end

%%%
%
function f = function_f(x1,x2)
%
f = (1-x1)^2 + 100*(x2 - x1^2)^2;
return
end

