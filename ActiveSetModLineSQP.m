function [output1,output2] = ActiveSetModLineSQP(objective,x0,lam0,constraints,toler,itermax,neff,tau,roh)
%% Active Set Line Search SQP Algorithm
%With Penalty Method Line Search
%% Source:
%Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (Second ed.). Springer. 
%doi:10.1007/978-0-387-40065-5
%Chapter 16: Quadratic Programming - Active Set Method for Convex QP's
%Chapter 18: Sequential Quadratic Programming - Intro and Merit Functions

%% Instructions for use:
%Define all functions as anonymous functions:
%funck - your objective function
%c(#)k - your constraint functions
%CHANGE lam0 to initial guess for lagrange multipliers (column vector)
%CHANGE size of lam0 to number of initial constraints
%CHANGE x0 to initial guess (column vector)
%To make a separate, non-function script, remove code from function definition, place
%into new script, and initialize the inputs listed below.

%% INPUTS
%objective: Objective Function, input as an anonymous function @(x) f(x)
%x0: Initial guess for states, input as column vector
%lam0: Initial guess for lagrange multipliers, input as column vector  
%constraints: Input as a vector function of anonymous functions @(x) [c1(x); c2(x);... cN(x)] 
%toler: Tolerance, typically set between 10E-4 and 10E-6
%itermax: Maximum number of iterations, typically set between 30 and 100
%neff: A constant parameter used in a penalty induced line search, set
%between 0.45 and 0.55
%tau: A constant parameter used to decrease step size, (0, 1), typically
%0.95
%roh: A constant parameter used to calculate the penalty, mu, typically set beween
%0.85 and 0.95

%% OUTPUTS
%output1: The state vector corresponding to a local or global minimum
%output2: A matrix of logicals showing which constraints are active at each
%iteration

%% Initialize SQP and Functions
%function f(x)
tic
funck = objective; %(for testing place an objective function here as shown) @(x) f(x(1),x(2))
%constraints in order 
%(for testing place constraints here as shown)
%c1k = @(x) x(1)-2*x(2)+2;
%c2k = @(x) -x(1)-2*x(2)+6;
%c3k = @(x) -x(1)+2*x(2)+2;
%c4k = @(x) x(1);
%c5k = @(x) x(2);
%function of constraints 
Cfunc = constraints; %(for testing create function as shown) @(x) [c1k(x); c2k(x); c3k(x); c4k(x)]; %c4k(x); c5k(x)];

%% Initial Guess: x0 and lagrange multipliers
x0 = x0;
lam0 = lam0;

%% Parameter Initialization
numc = length(Cfunc(x0)); %number of constraints
v = length(x0); %number of variables;
tol = toler;%10E-8;
constol = 10E-6; %constraint tolerance;
C = zeros(numc,1);
J = zeros(numc,v); %size of constraint jacobian
e = 0.05; %tolerance for finite difference methods
E = eye(v);
cutval = 1E-9;
itermax = itermax;%1000;
itermain = 1;
neff = neff;%0.45;
tau = tau;%0.95;
roh = roh;%0.95;
sigma = 1;
nullactset = 0;
initialflag = 0;
blockflag = 0;

%% Linearization Step
k = 1;
while itermain <= itermax
m = 1;
n = 1;
%create jacobian of constraints
while n <= v
    while m <= numc
        index = @(Cfunc,i) Cfunc(i);
        J(m,n) = (index(Cfunc((x0'+e*E(n,:))),m)-index(Cfunc((x0')),m))./e;
        m = m+1;
    end
    m = 1;
    n = n+1;
end

%% Alter Constraints Vector
count = 1;
while count <= numc
    index = @(Cfunc,i) Cfunc(i);
    C(count) = index(Cfunc(x0),count);
    count = count+1;
end

%% Finite Difference Lagrangian
Lx = @(x,lam0) funck(x)-lam0'*[Cfunc(x)]; %c3k(x)];
%Initialize Hessian of Lagrangian as I(v);
del2L = eye(v);
m = 1;
n = 1;
while m <= v
    while n <= v
        del2L(m,n) = (Lx(x0'+e*E(m,:)+e*E(n,:),lam0)-Lx(x0'+e*E(m,:),lam0)-Lx(x0'+e*E(n,:),lam0)+Lx(x0',lam0))./(e^2);
        if del2L(m,n) <= cutval
            del2L(m,n) = 0;
        end
        n = n+1;
    end
    n = 1;
    m = m+1;
end
fk = funck(x0'); %function value of f at x
n = 1;
delfk = zeros(v,1); %gradient of f at x, equivalence to gradient of Lagrangian holds due to constraints A(x)*p+c(x) = 0;
while n <= v
    delfk(n) = (funck(x0'+e*E(n,:))-funck(x0'-e*E(n,:)))./(2*e);
    n = n+1;
end

%% Initialize QP
x_0 = x0;
Bk = del2L;
A = J;
b = -C;

%% Find Active Set and Form KKT

if initialflag == 0
    initialflag = 1;
    ActConst(:,k) = zeros(length(b),1);
    constraintval = Cfunc(x_0); %Warm Starts Solution
    citer = 1;
    while citer <= length(constraintval)
        if constraintval(citer) > constol
            ActConst(citer,k) = 0;
        end
        if constraintval(citer) <= constol
            ActConst(citer,k) = 1;
        end
        citer = citer+1;
    end
    if sum(ActConst(:,k)) == 0
        nullactset = 1;
    else
        WorksetA = A(logical(ActConst(:,k)),:);
        Worksetb = b(logical(ActConst(:,k)),:);
    end
%{
else %Modify the Active Set Logical Placeholder for Subsequent Iterations
    if nullactset ~= 1 %change zerocheck to nullactset
        [~,modconst] = ismember(Worksetb,round(b));
        ActConst((modconst),k) = 1;
    else
        ActConst(:,k) = zeros(length(b),1);
    end
%}
end
c = del2L*x0-delfk;
%% BEGIN QP LOOP
while k >= 1
if nullactset == 0
    if initialflag ~= 0
        Cnew = -Cfunc(x_0);
        Worksetb = Cnew(logical(ActConst(:,k)));
    end
    constrain = WorksetA*x_0-Worksetb;
    KKT = [Bk -WorksetA'; WorksetA zeros(length(Worksetb(:,1)))];
    cost = [c; constrain];
    negflag = 1;
else
    KKT = Bk;
    cost = c;
    negflag = 1;
end

xks = KKT\cost;
xx = xks(1:v);
pk = xx-x_0;

if length(xks) > v
    lamk = xks(v+1:length(xks));
else
    lamk = 0;
end
k = k+1;
%% Check for Descent Direction and Check Sign of Lagrange Multipliers 

if norm(pk) > -cutval && norm(pk) < cutval
    f = 1;
    i = 1;
    while i <= length(ActConst(:,1))
        if ActConst(i,k-1) == 1
            lam0(i) = lamk(f);
            f = f+1;
        end
        i = i+1;
    end
    neg = find(lamk(1:length(lamk))<0);
    sizeneg = numel(neg);
    if sizeneg > 0
        if length(neg) > 1
        testset = zeros(length(lamk),1);
        iter = 1;
        while iter <= length(lamk)
            if iter == neg(iter)
            testset(iter) = lamk(iter);
            else
            testset(iter) = 0;
            end
        iter = iter + 1;
        end
        [~,remove] = min(testset);
        index = true(length(WorksetA(:,1)), size(length(WorksetA(:,1)),1));
        index(remove) = false;
        % Remove Assosciated Constraint from Working Set
        WorksetA = WorksetA(index,:);
        Worksetb = Worksetb(index,:);
        newact = find(lam0 == min(testset));
        ActConst(:,k) = ActConst(:,k-1);
        ActConst(newact,k) = 0;
        else
            zerocheck = 1;
            nullactset = 1;
            ActConst(:,k) = zeros(length(ActConst(:,k-1)),1);
        end
    else
        x_0 = x_0;
        break
    end
else
    i = 1;
    f = 1;
    lamnew = zeros(length(lam0),1);
while i <= length(ActConst(:,k-1))
    if ActConst(i,k-1) == 1
        lamnew(i) = lamk(f);
        f = f+1;
    else
        lamnew(i) = lam0(i);
    end
    i = i+1;
end
lamk = lamnew;
plam = lamk-lam0;
blockcheck = -Cfunc(x_0)./(A*pk);
biter = 1;
while biter <= length(blockcheck)
    if blockcheck(biter) > tol && blockcheck(biter) <= 1
        valid(biter) = blockcheck(biter);
    else
        valid(biter) = 2;
    end
    biter = biter + 1;
end
newalpha = min(valid);

if newalpha == 1
    %nullactset = 1; %?
    blockflag = 0;
    ActConst(:,k) = ActConst(:,k-1);
    x_0 = x_0+newalpha*pk;
end
if newalpha > 0 && newalpha < 1
   block = find(blockcheck == newalpha);
   WorksetA = A(block,:);
   Worksetb = b(block,:);
   ActConst(:,k) = zeros(length(blockcheck),1);
   ActConst(block,k) = 1;
   x_0 = x_0+newalpha*pk;
   nullactset = 0;
end
if newalpha > 1
   x_0 = x_0+pk;
   ActConst(:,k) = ActConst(:,k-1);
end
end
end
%% END QP LOOP

%% Determine Step Size Using Modified Line Search, Enforce Convergence
pkstar = x_0-x0;
plamstar = lamnew-lam0;
Ck = Cfunc(x0);
mu = (delfk'*pkstar+(sigma/2)*pkstar'*del2L*pkstar)/((1-roh)*norm(Ck,1));%(delfk'*pk+(sigma/2)*pk'*del2L*pk)/((1-roh)*norm(C));
alpha = 1;
while (funck(x0+alpha*pkstar)+mu*norm(Cfunc(x0+alpha*pkstar),1)) > (funck(x0)+mu*norm(Cfunc(x0))+neff*alpha*(delfk'*pkstar-mu*norm(Ck,1)))%(funck(x_0)+mu*norm(Cfunc(x_0))+neff*alpha*(-pk'*del2L*pk-(mu-norm(lamnew))*norm(Cfunc(x_0))))%Lx(x_0+alpha*pk,lam0+alpha*plam)+mu*norm(Cfunc(x_0+alpha*pk)) > (Lx(x_0,lam0)+mu*norm(Cfunc(x_0))+neff*alpha*(-pk'*del2L*pk-(mu-norm(lamnew))*norm(Cfunc(x_0))))
        alpha = tau*alpha;
end
xk = x0+alpha*pkstar;
lam0k = lam0+alpha*plamstar;
itermain = itermain + 1;

%% Tolerance Test
if abs(funck(xk)-funck(x0)) <= tol
    break
end
x0 = x0+alpha*pkstar;
lam0 = lam0+alpha*plamstar;
initialflag = 0;
end
output1 = x0;
output2 = ActConst;
toc

end