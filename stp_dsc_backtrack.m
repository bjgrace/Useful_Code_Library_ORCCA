function [ak] = stp_dsc_backtrack(aini,x1,p,c)
%% Backtracking Line Search
%A simple backtracking line search obeying the Wolfe conditions

%% Source
%Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (Second ed.). Springer. 
%doi:10.1007/978-0-387-40065-5

%% Inputs
%aini: Initial guess for a
%x1: Current value of x
%p: Constant to reduce value of alpha (0,1);
%c: Constant to lead to sufficient decrease (0,1];
%define objective function f here

%% Outputs
%ak: New step size

%% Function Definition
f = @(x) -100*exp(-0.5*x)*sin(-2*x);
%gradient
g = @(x) -100*exp(-0.5*x)*(-0.5*sin(-2*x)-2*cos(-2*x));
%1D Hessian (2nd derivative)
B = @(x) 50*exp(-0.5*x)*((15/2)*sin(-2*x)-4*cos(-2*x));

%% Evaluate Function, Derivative, Hessian at x_k
fk = f(x1); %equivalent to phi(0)
gk = g(x1); %equivalent to phi'(0)
Bk = B(x1);

%% Find step direction
pk = -(1/Bk)*gk;
a = aini;
while f(x1 + a*pk) <= fk+(c*a*gk*pk)
    a = p*a;
    fk = f(x1 + a*pk);
    gk = g(x1 + a*pk);
    if a <= 0
        break
    end
end
ak = a;
end




