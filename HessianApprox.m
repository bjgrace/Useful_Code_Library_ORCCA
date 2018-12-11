function [hessian,symhessian] = HessianApprox(nv,x0,coefficients,const)
%% Note: For QP Problem Formulation
%% INPUTS
%coefficients: Array of coefficients of QP in order - i.e. if problem is
%C1*x1^2+C2*x2^2+F, then coeff is [C1 0 C2 0];
%const: Array of constants - i.e. if problem is above, const = [F];
%x0: Current x value
%nv: Number of variables, set to length(x0)

%% OUTPUTS
%hessian: The Hessian Matrix
%symhessian: A reformulated, symmetric Hessian Matrix

coeff = coefficients;
co = sum(const);
v = nv;
x_0 = x0;
e = 0.5;
i = 1;
m = 1;
n = 1;
A = zeros(v);
C = coeff((v.^2+1):length(coeff));
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
while i <= v.^2
    if n > v
        n = 1;
        m = m + 1;
    end
    B(m,n) = (x_0+e*E(:,m)+e*E(:,n))'*A*(x_0+e*E(:,m)+e*E(:,n))+C*(x_0+e*E(:,m)+e*E(:,n))+co;
    B(m,n) = B(m,n)-((x_0+e*E(:,m))'*A*(x_0+e*E(:,m))+C*(x_0+e*E(:,m))+co);
    B(m,n) = B(m,n)-((x_0+e*E(:,n))'*A*(x_0+e*E(:,n))+C*(x_0+e*E(:,n))+co);
    B(m,n) = B(m,n)+((x_0)'*A*(x_0)+C*(x_0)+co);
    B(m,n) = B(m,n)/(e.^2);
    i = i + 1;
    n = n + 1;
end
hessian = B;
symhessian = B'*B;