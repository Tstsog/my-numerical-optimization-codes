% This Matlab code performs an optimization (finding the maximum/minimum) for one-variable function f(x) 
% using Newton-Raphson iterative method [1].
%
% Ref. [1] S. Chapra, "Applied numerical methods with MATLAB", Mc Craw Hill, Singapore (2008).
%
% The non-variable function: f(x) = 2*sin(x) - x^2.  
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
function [] = one_variable_minimization_newton_1
clc; clear one_variable_minimization_newton_1; 
%
format long
%
x = 0.500; tol = 1e-6;
dx = 0.0001; iter_max = 50.;
%
for it = 1:iter_max
    %
    f = functon_1_var(x);
    %
    output = [it, x, f]
    % central difference scheme
    Jacobian = (functon_1_var(x+dx) - functon_1_var(x-dx))/(2*dx);
    Hessian = (functon_1_var(x+dx) - 2*functon_1_var(x) + functon_1_var(x-dx))/dx^2;        
    %
    % forward difference scheme 
%    Jacobian = (functon_1_var(x+dx) - functon_1_var(x))/dx;    
%    Hessian = (functon_1_var(x+2*dx) - 2*functon_1_var(x+dx) + functon_1_var(x))/dx^2;        
    %
    x = x - Hessian\Jacobian;
    %
    if (abs(Jacobian) <= tol)
        break; 
    end
    %
end
%
[Jacobian, Hessian]; % -0.000000000016265  -3.347223837324442 < 0 (local maximum)
%
% output = [it, x, f]
% output = 1.000000000000000   2.000000000000000  -2.181405146348637
% output = 2.000000000000000   0.734536155144740   0.800942615050591
% output = 3.000000000000000   0.739089650329804   0.800977224192609
% output = 4.000000000000000   0.739085059618187   0.800977224226749


%%%
figure(1)
hold on
plot(output(2), output(3), 'r+', 'MarkerSize', 12, LineWidth=1.5 )
plot((0.:0.01:1.5), 2.*sin((0.:0.01:1.5)) - (0.:0.01:1.5).^2, 'b', LineWidth=1.5 )
hold off
xlabel('$x$','interpreter','latex')
ylabel('$f(x)$','interpreter','latex', 'Rotation', 1)
set(gca,'FontSize',16)
box on

%%%
return
end
%
function f = functon_1_var(x)
%
f = 2*sin(x) - x.^2;
%%%
return
end
