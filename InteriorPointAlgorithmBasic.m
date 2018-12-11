%% Interior Point Algorithm for Bounded, Equality, and Inequality Constrained Problems
%% Notes:
%{
-Step size calculation could use refinement, current perturbation to step
seems to work, code for more extensive calculation is commented out and
included in step size calculation section
-Unlike SQP method, does not use the estimated control inputs as an initial guess,
just the state so that inputs are not out of bounds of inequality
constrained region
%}
%% Source:
%Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (Second ed.). Springer. 
%doi:10.1007/978-0-387-40065-5
%Chapter 19: Interior Point Methods (non-linear)

%% INPUTS/Changeable Variables
%{
objmat: The block diagonal matrix of weights, Q, R, and S, used in the objective
function
tol: Tolerance, set between 10E-4 and 10E-6
maxiter: Maximum number of iterations
N: Horizon length
v: length of state vector
m: length of input vector
delT: Timestep
st: Current combined state/input vector
x0: Initial state (current state after 1st iteration of model) 
u0: Initial input
I: Moment of inertia matrix
simga: Constant parameter (0,1) typically between 0.5 and 0.95
mu: Barrier parameter mu > 0
umin: Minimum allowable input
umax: Maximum allowable input
%}

%% OUTPUTS
%{
L: The control input torque calculated for the 1st move along
the time horizon
%}
% function definition at later point in time
%% Problem Initialization
I = [0.0833 0 0; 0 0.1083 0; 0 0 0.0417];
I1 = I(1,1);
I2 = I(2,2);
I3 = I(3,3);
delT = 0.1;
N = 5;
v = 6;
m = 3;
tol = 10E-4;
sigma = 0.95;
maxiter = 30;
yflag = 1;
zflag = 1;
bfgsflag = 0;
mu = 2;
tau = 0.995;
x0 = x0;%ones((N+1)*v,1);
s0 = ones(2*N*m,1);
s = s0;
u0 = zeros(N*m,1);
lam0 = ones(((N+1)*v),1);
lam20 = ones((2*N*m),1);
z0 = lam20;
y = lam0;
z = lam20;
S = zeros(length(s));
Z = zeros(length(z));
e = ones(length(s),1);
umin = -0.05;
umax = 0.05;
objmat = objmat;
x = [x0; zeros((N)*v+N*m,1)];
%x = [x0;u0];
% Preallocation
Ceqfunc = zeros(((N+1)*v),1);
Cineqfunc = zeros((2*N*m),1);
Jeq = zeros((N+1)*v,(N+1)*v+N*m);
Jineq = zeros(2*N*m,(N+1)*v+N*m);
Heq = zeros(length(x));
Hineq = zeros(length(x));
xf = zeros(length(x),1);
xb = zeros(length(x),1);
xf1 = zeros(length(x),1);
xb1 = zeros(length(x),1);
xb2 = zeros(length(x),1);
%% Assign slack variables to S
i = 1;
j = 1;
while i <= length(s)
    while j <= length(s)
        if i == j
            S(i,j) = s(i);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end
%% Calculate z0 (first iter only) and assign slack variables
if zflag == 1
    z = S\(mu*e);
    zflag = 0;
end
i = 1;
j = 1;
while i <= length(z)
    while j <= length(z)
        if i == j
            Z(i,j) = z(i);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end
%% Propogate Un-perturbed/Un-controlled Dynamics
propstore = zeros(N*v,1);
i = 1;
while i <= N
    if i == 1
        propstore(1:3) = (1/2)*[x0(1) sqrt(abs(1-x0(1).^2-x0(2).^2-x0(3).^2)) -x0(3) x0(2); x0(2) x0(3) sqrt(abs(1-x0(1).^2-x0(2).^2-x0(3).^2)) -x0(1); x0(3) -x0(2) x0(1) sqrt(abs(1-x0(1).^2-x0(2).^2-x0(3).^2))]*[0; x0(4); x0(5); x0(6)]*delT;
        propstore(4:6) = inv(I)*(-[0 -x0(6) x0(5); x0(6) 0 -x0(4); -x0(5) x0(4) 0]*I*[x0(4); x0(5); x0(6)]*delT);
    else
        propstore((i-1)*v+1:(i-1)*v+3) = (1/2)*[propstore((i-1)*1) sqrt(abs(1-propstore((i-1)*1).^2-propstore((i-1)*2).^2-propstore((i-1)*3).^2)) -propstore((i-1)*3) propstore((i-1)*2); propstore((i-1)*2) propstore((i-1)*3) sqrt(abs(1-propstore((i-1)*1).^2-propstore((i-1)*2).^2-propstore((i-1)*3).^2)) -propstore((i-1)*1); propstore((i-1)*3) -propstore((i-1)*2) propstore((i-1)*1) sqrt(abs(1-propstore((i-1)*1).^2-propstore((i-1)*2).^2-propstore((i-1)*3).^2))]*[0; propstore((i-1)*4); propstore((i-1)*5); propstore((i-1)*6)]*delT;
        propstore((i-1)*v+4:(i-1)*v+6) = inv(I)*(-[0 -propstore((i-1)*6) propstore((i-1)*5); propstore((i-1)*6) 0 -propstore((i-1)*4); -propstore((i-1)*5) propstore((i-1)*4) 0]*I*[propstore((i-1)*4); propstore((i-1)*5); propstore((i-1)*6)])*delT;
    end
    i = i+1;
end
x(v+1:(N+1)*v) = propstore;
iter = 1;
true = 1;
while iter <= maxiter
while true >= 1
%% Equality Constraint Function Values/Inequality Constraint Function Values
Ceqfunc(1:6) = [x(1)-x(1); x(2)-x(2); x(3)-x(3); x(4)-x(4); x(5)-x(5); x(6)-x(6)];
i = 1;
while i <= N
    Ceqfunc(7+(i-1)*v:6+i*v) = [x(7+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))-x(3+(i-1)*v)*x(5+(i-1)*v)+x(2+(i-1)*v)*x(6+(i-1)*v))*delT;
                      x(8+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))+x(3+(i-1)*v)*x(4+(i-1)*v)-x(1+(i-1)*v)*x(6+(i-1)*v))*delT;
                      x(9+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))+x(1+(i-1)*v)*x(5+(i-1)*v)-x(2+(i-1)*v)*x(4+(i-1)*v))*delT;
                      x(10+(i-1)*v)-((I2*x(5+(i-1)*v)*x(6+(i-1)*v)-I3*x(5+(i-1)*v)*x(6+(i-1)*v)+x(25+(i-1)*m))/I1)*delT;
                      x(11+(i-1)*v)-((-I1*x(4+(i-1)*v)*x(6+(i-1)*v)+I3*x(4+(i-1)*v)*x(6+(i-1)*v)+x(26+(i-1)*m))/I2)*delT;
                      x(12+(i-1)*v)-((I1*x(4+(i-1)*v)*x(5+(i-1)*v)-I2*x(4+(i-1)*v)*x(5+(i-1)*v)+x(27+(i-1)*m))/I3)*delT;];
                  i = i+1;
end
i = 1;
while i <= N
    Cineqfunc(1+(i-1)*2*m:i*(2*m)) = [x((N+1)*v+1+(i-1)*m)-umin; x((N+1)*v+2+(i-1)*m)-umin; x((N+1)*v+3+(i-1)*m)-umin;
                                      -x((N+1)*v+1+(i-1)*m)+umax; -x((N+1)*v+2+(i-1)*m)+umax; -x((N+1)*v+3+(i-1)*m)+umax;];%[umin; umin; umin; -umax; -umax; -umax];
    i = i+1;
end
%% Equality Constraint Jacobian Values/Inequality Constraint Jacobian Values
E = eye(length(x));
efd = 0.05;
forw = zeros(N*v,1);
backw = zeros(N*v,1);
forwineq = zeros(N*2*m,1);
backineq = zeros(N*2*m,1);
i = 1;
j = 1;
while i <= N
    while j <= length(x)
    xf = x+efd*E(:,j);
    xb = x-efd*E(:,j);
    forw = [xf(7+(i-1)*v)-(1/2)*(sqrt(abs(1-xf(1+(i-1)*v)^2-xf(2+(i-1)*v)^2-xf(3+(i-1)*v)^2))-xf(3+(i-1)*v)*xf(5+(i-1)*v)+xf(2+(i-1)*v)*xf(6+(i-1)*v))*delT;
                      xf(8+(i-1)*v)-(1/2)*(sqrt(abs(1-xf(1+(i-1)*v)^2-xf(2+(i-1)*v)^2-xf(3+(i-1)*v)^2))+xf(3+(i-1)*v)*xf(4+(i-1)*v)-xf(1+(i-1)*v)*xf(6+(i-1)*v))*delT;
                      xf(9+(i-1)*v)-(1/2)*(sqrt(abs(1-xf(1+(i-1)*v)^2-xf(2+(i-1)*v)^2-xf(3+(i-1)*v)^2))+xf(1+(i-1)*v)*xf(5+(i-1)*v)-xf(2+(i-1)*v)*xf(4+(i-1)*v))*delT;
                      xf(10+(i-1)*v)-((I2*xf(5+(i-1)*v)*xf(6+(i-1)*v)-I3*xf(5+(i-1)*v)*xf(6+(i-1)*v)+xf(25+(i-1)*m))/I1)*delT;
                      xf(11+(i-1)*v)-((-I1*xf(4+(i-1)*v)*xf(6+(i-1)*v)+I3*xf(4+(i-1)*v)*xf(6+(i-1)*v)+xf(26+(i-1)*m))/I2)*delT;
                      xf(12+(i-1)*v)-((I1*xf(4+(i-1)*v)*xf(5+(i-1)*v)-I2*xf(4+(i-1)*v)*xf(5+(i-1)*v)+xf(27+(i-1)*m))/I3)*delT;];
    backw = [xb(7+(i-1)*v)-(1/2)*(sqrt(abs(1-xb(1+(i-1)*v)^2-xb(2+(i-1)*v)^2-xb(3+(i-1)*v)^2))-xb(3+(i-1)*v)*xb(5+(i-1)*v)+xb(2+(i-1)*v)*xb(6+(i-1)*v))*delT;
                      xb(8+(i-1)*v)-(1/2)*(sqrt(abs(1-xb(1+(i-1)*v)^2-xb(2+(i-1)*v)^2-xb(3+(i-1)*v)^2))+xb(3+(i-1)*v)*xb(4+(i-1)*v)-xb(1+(i-1)*v)*xb(6+(i-1)*v))*delT;
                      xb(9+(i-1)*v)-(1/2)*(sqrt(abs(1-xb(1+(i-1)*v)^2-xb(2+(i-1)*v)^2-xb(3+(i-1)*v)^2))+xb(1+(i-1)*v)*xb(5+(i-1)*v)-xb(2+(i-1)*v)*xb(4+(i-1)*v))*delT;
                      xb(10+(i-1)*v)-((I2*xb(5+(i-1)*v)*xb(6+(i-1)*v)-I3*xb(5+(i-1)*v)*xb(6+(i-1)*v)+xb(25+(i-1)*m))/I1)*delT;
                      xb(11+(i-1)*v)-((-I1*xb(4+(i-1)*v)*xb(6+(i-1)*v)+I3*xb(4+(i-1)*v)*xb(6+(i-1)*v)+xb(26+(i-1)*m))/I2)*delT;
                      xb(12+(i-1)*v)-((I1*xb(4+(i-1)*v)*xb(5+(i-1)*v)-I2*xb(4+(i-1)*v)*xb(5+(i-1)*v)+xb(27+(i-1)*m))/I3)*delT;];
    forwineq = [xf((N+1)*v+1+(i-1)*m)-umin; xf((N+1)*v+2+(i-1)*m)-umin; xf((N+1)*v+3+(i-1)*m)-umin; -xf((N+1)*v+1+(i-1)*m)+umax; -xf((N+1)*v+2+(i-1)*m)+umax; -xf((N+1)*v+3+(i-1)*m)+umax;];
    backineq = [xb((N+1)*v+1+(i-1)*m)-umin; xb((N+1)*v+2+(i-1)*m)-umin; xb((N+1)*v+3+(i-1)*m)-umin; -xb((N+1)*v+1+(i-1)*m)+umax; -xb((N+1)*v+2+(i-1)*m)+umax; -xb((N+1)*v+3+(i-1)*m)+umax;];
    Jeq(7+(i-1)*v:6+i*v,j) = (forw-backw)./(2*efd);
    Jineq(1+(i-1)*2*m:i*2*m,j) = (forwineq-backineq)./(2*efd);
    j = j+1;
    end
    j = 1;
    i = i+1;
end
Jeq(1:6,1:6) = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
%% Cost Function
objfunc = x'*objmat*x;
%% Cost Function Gradient
delf = zeros((N+1)*v+N*m,1);
i = 1;
while i <= length(x)
    xf = x+efd*E(:,i);
    delf(i) = xf'*objmat*xf;
    i = i+1;
end

%% Cost Function Hessian
Hfunc = zeros((N+1)*v+N*m);
i = 1;
j = 1;
while i <= length(x)
    while j <= length(x)
        xf1 = x+efd*E(:,i)+efd*E(:,j);
        xb1 = x+efd*E(:,i);
        xb2 = x+efd*E(:,j);
        Hfunc(i,j) = (xf1'*objmat*xf1+x'*objmat*x-(xb1'*objmat*xb1)-(xb2'*objmat*xb2))./(efd^2);
        j = j+1;
    end
    j = 1;
    i = i+1;
end
%% Calculate initial y (1st iteration only)
if yflag == 1
    y0 = (delf-Jineq'*z0)\Jeq';
    yflag = 0;
    y = y0';
end
%% Lagrangian Gradient
delL = delf-Jeq'*y-Jineq'*z;
%% Equality Constraint, Inequality Constraint Hessian Values
i = 1;
j = 1;
k = 1;
if bfgsflag == 0
forw1 = zeros((N+1)*v,1);
backw1 = zeros((N+1)*v,1);
backw2 = zeros((N+1)*v,1);
center = zeros((N+1)*v,1);
forwineq1 = zeros(2*N*m,1);
backineq1 = zeros(2*N*m,1);
backineq2 = zeros(2*N*m,1);
centerineq = zeros(2*N*m,1);
while j <= length(x)
     while k <= length(x)
         xf1 = x+efd*E(:,j)+efd*E(:,k);
         xb1 = x+efd*E(:,j);
         xb2 = x+efd*E(:,k);
         while i <= N
         %f1 + f0 - fb1 - fb2 / e^2
            forw1(7+(i-1)*v:6+i*v) = [xf1(7+(i-1)*v)-(1/2)*(sqrt(abs(1-xf1(1+(i-1)*v)^2-xf1(2+(i-1)*v)^2-xf1(3+(i-1)*v)^2))-xf1(3+(i-1)*v)*xf1(5+(i-1)*v)+xf1(2+(i-1)*v)*xf1(6+(i-1)*v))*delT;
                   xf1(8+(i-1)*v)-(1/2)*(sqrt(abs(1-xf1(1+(i-1)*v)^2-xf1(2+(i-1)*v)^2-xf1(3+(i-1)*v)^2))+xf1(3+(i-1)*v)*xf1(4+(i-1)*v)-xf1(1+(i-1)*v)*xf1(6+(i-1)*v))*delT;
                   xf1(9+(i-1)*v)-(1/2)*(sqrt(abs(1-xf1(1+(i-1)*v)^2-xf1(2+(i-1)*v)^2-xf1(3+(i-1)*v)^2))+xf1(1+(i-1)*v)*xf1(5+(i-1)*v)-xf1(2+(i-1)*v)*xf1(4+(i-1)*v))*delT;
                   xf1(10+(i-1)*v)-((I2*xf1(5+(i-1)*v)*xf1(6+(i-1)*v)-I3*xf1(5+(i-1)*v)*xf1(6+(i-1)*v)+xf1(25+(i-1)*m))/I1)*delT;
                   xf1(11+(i-1)*v)-((-I1*xf1(4+(i-1)*v)*xf1(6+(i-1)*v)+I3*xf1(4+(i-1)*v)*xf1(6+(i-1)*v)+xf1(26+(i-1)*m))/I2)*delT;
                   xf1(12+(i-1)*v)-((I1*xf1(4+(i-1)*v)*xf1(5+(i-1)*v)-I2*xf1(4+(i-1)*v)*xf1(5+(i-1)*v)+xf1(27+(i-1)*m))/I3)*delT;];
            backw1(7+(i-1)*v:6+i*v) = [xb1(7+(i-1)*v)-(1/2)*(sqrt(abs(1-xb1(1+(i-1)*v)^2-xb1(2+(i-1)*v)^2-xb1(3+(i-1)*v)^2))-xb1(3+(i-1)*v)*xb1(5+(i-1)*v)+xb1(2+(i-1)*v)*xb1(6+(i-1)*v))*delT;
                   xb1(8+(i-1)*v)-(1/2)*(sqrt(abs(1-xb1(1+(i-1)*v)^2-xb1(2+(i-1)*v)^2-xb1(3+(i-1)*v)^2))+xb1(3+(i-1)*v)*xb1(4+(i-1)*v)-xb1(1+(i-1)*v)*xb1(6+(i-1)*v))*delT;
                   xb1(9+(i-1)*v)-(1/2)*(sqrt(abs(1-xb1(1+(i-1)*v)^2-xb1(2+(i-1)*v)^2-xb1(3+(i-1)*v)^2))+xb1(1+(i-1)*v)*xb1(5+(i-1)*v)-xb1(2+(i-1)*v)*xb1(4+(i-1)*v))*delT;
                   xb1(10+(i-1)*v)-((I2*xb1(5+(i-1)*v)*xb1(6+(i-1)*v)-I3*xb1(5+(i-1)*v)*xb1(6+(i-1)*v)+xb1(25+(i-1)*m))/I1)*delT;
                   xb1(11+(i-1)*v)-((-I1*xb1(4+(i-1)*v)*xb1(6+(i-1)*v)+I3*xb1(4+(i-1)*v)*xb1(6+(i-1)*v)+xb1(26+(i-1)*m))/I2)*delT;
                   xb1(12+(i-1)*v)-((I1*xb1(4+(i-1)*v)*xb1(5+(i-1)*v)-I2*xb1(4+(i-1)*v)*xb1(5+(i-1)*v)+xb1(27+(i-1)*m))/I3)*delT;];
            backw2(7+(i-1)*v:6+i*v) = [xb2(7+(i-1)*v)-(1/2)*(sqrt(abs(1-xb2(1+(i-1)*v)^2-xb2(2+(i-1)*v)^2-xb2(3+(i-1)*v)^2))-xb2(3+(i-1)*v)*xb2(5+(i-1)*v)+xb2(2+(i-1)*v)*xb2(6+(i-1)*v))*delT;
                   xb2(8+(i-1)*v)-(1/2)*(sqrt(abs(1-xb2(1+(i-1)*v)^2-xb2(2+(i-1)*v)^2-xb2(3+(i-1)*v)^2))+xb2(3+(i-1)*v)*xb2(4+(i-1)*v)-xb2(1+(i-1)*v)*xb2(6+(i-1)*v))*delT;
                   xb2(9+(i-1)*v)-(1/2)*(sqrt(abs(1-xb2(1+(i-1)*v)^2-xb2(2+(i-1)*v)^2-xb2(3+(i-1)*v)^2))+xb2(1+(i-1)*v)*xb2(5+(i-1)*v)-xb2(2+(i-1)*v)*xb2(4+(i-1)*v))*delT;
                   xb2(10+(i-1)*v)-((I2*xb2(5+(i-1)*v)*xb2(6+(i-1)*v)-I3*xb2(5+(i-1)*v)*xb2(6+(i-1)*v)+xb2(25+(i-1)*m))/I1)*delT;
                   xb2(11+(i-1)*v)-((-I1*xb2(4+(i-1)*v)*xb2(6+(i-1)*v)+I3*xb2(4+(i-1)*v)*xb2(6+(i-1)*v)+xb2(26+(i-1)*m))/I2)*delT;
                   xb2(12+(i-1)*v)-((I1*xb2(4+(i-1)*v)*xb2(5+(i-1)*v)-I2*xb2(4+(i-1)*v)*xb2(5+(i-1)*v)+xb2(27+(i-1)*m))/I3)*delT;];
            center(7+(i-1)*v:6+i*v) = [x(7+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))-x(3+(i-1)*v)*x(5+(i-1)*v)+x(2+(i-1)*v)*x(6+(i-1)*v))*delT;
                   x(8+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))+x(3+(i-1)*v)*x(4+(i-1)*v)-x(1+(i-1)*v)*x(6+(i-1)*v))*delT;
                   x(9+(i-1)*v)-(1/2)*(sqrt(abs(1-x(1+(i-1)*v)^2-x(2+(i-1)*v)^2-x(3+(i-1)*v)^2))+x(1+(i-1)*v)*x(5+(i-1)*v)-x(2+(i-1)*v)*x(4+(i-1)*v))*delT;
                   x(10+(i-1)*v)-((I2*x(5+(i-1)*v)*x(6+(i-1)*v)-I3*x(5+(i-1)*v)*x(6+(i-1)*v)+x(25+(i-1)*m))/I1)*delT;
                   x(11+(i-1)*v)-((-I1*x(4+(i-1)*v)*x(6+(i-1)*v)+I3*x(4+(i-1)*v)*x(6+(i-1)*v)+x(26+(i-1)*m))/I2)*delT;
                   x(12+(i-1)*v)-((I1*x(4+(i-1)*v)*x(5+(i-1)*v)-I2*x(4+(i-1)*v)*x(5+(i-1)*v)+x(27+(i-1)*m))/I3)*delT;];
            forwineq1(1+(i-1)*2*m:i*(2*m)) = [xf1((N+1)*v+1+(i-1)*m)-umin; xf1((N+1)*v+2+(i-1)*m)-umin; xf1((N+1)*v+3+(i-1)*m)-umin; -xf1((N+1)*v+1+(i-1)*m)+umax; -xf1((N+1)*v+2+(i-1)*m)+umax; -xf1((N+1)*v+3+(i-1)*m)+umax;];
            backineq1(1+(i-1)*2*m:i*(2*m)) = [xb1((N+1)*v+1+(i-1)*m)-umin; xb1((N+1)*v+2+(i-1)*m)-umin; xb1((N+1)*v+3+(i-1)*m)-umin; -xb1((N+1)*v+1+(i-1)*m)+umax; -xb1((N+1)*v+2+(i-1)*m)+umax; -xb1((N+1)*v+3+(i-1)*m)+umax;];
            backineq2(1+(i-1)*2*m:i*(2*m)) = [xb2((N+1)*v+1+(i-1)*m)-umin; xb2((N+1)*v+2+(i-1)*m)-umin; xb2((N+1)*v+3+(i-1)*m)-umin; -xb2((N+1)*v+1+(i-1)*m)+umax; -xb2((N+1)*v+2+(i-1)*m)+umax; -xb2((N+1)*v+3+(i-1)*m)+umax;];
            centerineq(1+(i-1)*2*m:i*(2*m)) = [x((N+1)*v+1+(i-1)*m)-umin; x((N+1)*v+2+(i-1)*m)-umin; x((N+1)*v+3+(i-1)*m)-umin; -x((N+1)*v+1+(i-1)*m)+umax; -x((N+1)*v+2+(i-1)*m)+umax; -x((N+1)*v+3+(i-1)*m)+umax;];
         i = i+1;
         end
         Heq(j,k) = (lam0'*(forw1+center-backw1-backw2))./(efd^2);
         Hineq(j,k) = (lam20'*(forwineq1+centerineq-backineq1-backineq2))./(efd^2);
         i = 1;
         k = k+1;
     end
     k = 1;
     j = j+1;
end
bfgsflag = 1;
%% Lagrangian Hessian
del2L = Hfunc-Heq-Hineq;
else
    sk = px;
    yk = delL-delLprev;
    del2L = del2L + (yk*yk')./(yk'*sk) - (del2L*sk*sk'*del2L)./(sk'*del2L*sk);
end
Error = max([norm(delL,1); norm(S*z-mu*e,1); norm(Ceqfunc,1); norm(Cineqfunc-s,1)]);
if Error <= mu
    break
end
delLprev = delL;
%% Form Newton Step
% fix zeros
KKT = [del2L zeros((N+1)*v+N*m,N*2*m) -Jeq' -Jineq'; zeros(2*N*m,(N+1)*v+N*m) Z zeros(N*2*m,(N+1)*v) S; Jeq zeros((N+1)*v,N*2*m) zeros((N+1)*v) zeros((N+1)*v,N*2*m); Jineq -eye(N*2*m) zeros(N*2*m,(N+1)*v) zeros(N*2*m)];
cost = -[delL; S*z-mu*e; Ceqfunc; Cineqfunc-s];

%% Newton Step
pvec = KKT\cost;
px = pvec(1:length(x));
ps = pvec(length(x)+1:length(x)+length(s));
py = pvec(length(x)+length(s)+1:length(x)+length(s)+length(y));
pz = pvec(length(x)+length(s)+length(y)+1:length(x)+length(s)+length(y)+length(z));

%% Calculate Step Size
alpha = 1;
%while norm(s+alpha*ps) <= norm((1-tau)*s)
%    alpha = alpha*tau;
%end
amaxs = 1;
amaxz = 1;
%{
amaxs = ((1-tau)*s-s)\ps;
if amaxs < 0
    amaxs = 0;
end
if amaxs > 1
    amaxs = 1;
end
%while norm(z+alpha*pz) <= norm((1-tau)*z)
%    alpha = alpha*tau;
%end
amaxz = ((1-tau)*z-z)\pz;
if amaxz < 0
    amaxz = 0;
end
if amaxz > 1
    amaxz = 1;
end
%}
%% Update Method
xOld = x;
sOld = s;
yOld = y;
zOld = z;
x = x+amaxs*px;
s = s+amaxs*ps;
y = y+amaxz*py;
z = z+amaxz*pz;
i = 1;
while i <= length(s)
    if s(i) <= 0
        s(i) = 0;
    end
    i = i+1;
end
i = 1;
while i <= length(z)
    if z(i) <= 0
        z(i) = 0;
    end
    i = i+1;
end
i = 1;
j = 1;
while i <= length(s)
    while j <= length(s)
        if i == j
            S(i,j) = s(i);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end
i = 1;
j = 1;
while i <= length(z)
    while j <= length(z)
        if i == j
            Z(i,j) = z(i);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end
true = true+1;
end
    if abs(x'*objmat*x-xOld'*objmat*xOld) <= tol
        break
    end
    mu = sigma*mu;
    iter = iter+1;
end
L = x((N+1)*v+1:(N+1)*v+3);
