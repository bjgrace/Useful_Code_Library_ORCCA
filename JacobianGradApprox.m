function [Jacobian,Gradient] = JacobianGradApprox(nv,x0,coefficients,const)
%% Note: For QP Problem Formulation
%% INPUTS
%coefficients: Array of coefficients of QP in order - i.e. if problem is
%C1*x1^2+C2*x2^2+F, then coeff is [C1 0 C2 0];
%const: Array of constants - i.e. if problem is above, const = [F];
%x0: Current x value
%nv: Number of variables, set to length(x0)

%% OUTPUTS
%Jacobian: The Jacobian
%Gradient: The Jacobian transpose, gradient of function

coeff = coefficients;
co = sum(const);
v = nv;
x_0 = x0;
e = 0.5;
i = 1;
m = 1;
n = 1;
A = zeros(v);
C = coeff(v.^2+1:v.^2+v);
while i <= v.^2
    if n > v
        n = 1;
        m = m + 1;
    end
    A(m,n) = coeff(i);
    i = i + 1;
    n = n + 1;
end
E = eye(v);
B = zeros(v);
i = 1;
m = 1;
n = 1;
while i <= v
    if n > v
        n = 1;
        m = m + 1;
    end
    J(m,n) = (x_0+e*E(:,n))'*A*(x_0+e*E(:,n))+C*(x_0+e*E(:,n))+co;
    J(m,n) = J(m,n)-((x_0)'*A*(x_0)+C*(x_0)+co);
    J(m,n) = J(m,n)./e;
    i = i + 1;
    n = n + 1;
end
Jacobian = J;
Gradient = Jacobian';