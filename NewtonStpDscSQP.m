%% SQP V2 Newton Step, Steepest Descent
%% Notes:
%{
-Currently not in function configuration
-Linearizes an objective function and constraints and forms QP for SQP and
solves
-Expanded upon in Simulink model to function as full SQP with function
definitions
%}

%% Instructions for use:
%Define all functions as anonymous functions:
%funck - your objective function
%c(#)k - your constraint functions
%CHANGE lam0 to initial guess for lagrange multipliers (1 is fine)
%CHANGE size of lam0 to number of constraints
%CHANGE numc to number of constraints
%CHANGE x0 to initial guess (column vector)
%Within Lx, add proper number of constraints (c(#)k)
%Within the alter constraints vector section, add proper number
%of C(#) = c(#)k
%%
%function f(x)
funck = @(x) 3*x(1).^2+2*x(1)*x(2)+x(1)*x(3)+2.5*x(2).^2+2*x(2)*x(3)+2*x(3).^2-8*x(1)-3*x(2)-3*x(3);%exp(x(1)*x(2)*x(3)*x(4)*x(5))-(1/2)*(x(1).^3 + x(2).^3 + 1).^2;
%constraints in order
c1k = @(x) x(1)+x(3)-3;%x(1).^2+x(2).^2+x(3).^2+x(4).^2+x(5).^2-10;
c2k = @(x) x(2)+x(3);%x(2)*x(3)-5*x(4)*x(5);
%c3k = @(x) x(1).^3+x(2).^3+1;

x0 = [0; 0; 0];%[-1.71; 1.59; 1.82; -0.763; -0.763;];
lam0 = [1; 1;]; %initial lagrange multiplier guess
numc = 2; %number of constraints
v = length(x0); %number of variables;

%initialization
tol = 10E-2;
C = zeros(numc,1);
A = zeros(numc,v); %size of constraint jacobian
e = 0.05; %tolerance for finite difference methods
E = eye(v);

itermax = 1000;
iter = 1;

while iter <= itermax
m = 1;
n = 1;
%create jacobian of constraints
while n <= v
    A(1,n) = (c1k(x0'+e*E(n,:))-c1k(x0'))./e;
    A(2,n) = (c2k(x0'+e*E(n,:))-c2k(x0'))./e;
    %A(3,n) = (c3k(x0'+e*E(n,:))-c3k(x0'))./e;
    n = n+1;
end
%% alter constraints vector
C(1) = c1k(x0);
C(2) = c2k(x0);
%C(3) = c3k(x0);
%%
%define Lagriangian as Lx
Lx = @(x,lam0) funck(x)-lam0'*[c1k(x); c2k(x);]; %c3k(x)];

%Initialize Hessian of Lagrangian as I(v);
del2L = eye(v);
m = 1;
n = 1;
while m <= v
    while n <= v
        del2L(m,n) = (Lx(x0'+e*E(m,:)+e*E(n,:),lam0)-Lx(x0'+e*E(m,:),lam0)-Lx(x0'+e*E(n,:),lam0)+Lx(x0',lam0))./(e^2);
        n = n+1;
    end
    n = 1;
    m = m+1;
end
Bk = del2L;

fk = funck(x0'); %function value of f at x

n = 1;
delfk = zeros(v,1); %gradient of f at x, equivalence to gradient of Lagrangian holds due to constraints A(x)*p+c(x) = 0;
while n <= v
    delfk(n) = (funck(x0'+e*E(n,:))-funck(x0'))./e;
    n = n+1;
end

KKT = [Bk -A'; A zeros(numc)];
cost = [-delfk; -C];
plam = KKT\cost;
pk = plam(1:v);
lam0 = plam(v+1:length(plam));
k = 1;
control = 0.8;
tau = 0.8;
alpha = (delfk'*delfk)/(delfk'*del2L*delfk);%0.0013;

xk = x0+alpha*pk;

x0 = xk;
xstore(iter,:) = x0;

iter = iter+1;
end
