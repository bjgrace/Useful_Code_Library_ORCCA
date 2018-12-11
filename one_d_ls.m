%% 1D Line Search (Basic Format, Wolfe Conditions, Curvature Conditions)
%Note to User: This function will result in an error if conditions are not
%satisfied, as alpha* will not be assigned a value. This script is highly
%dependent on initial guess for x as well as low and high bounds on the
%step size

%% Source
%Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (Second ed.). Springer. 
%doi:10.1007/978-0-387-40065-5

%% INPUTS
%a1: Low bound for step (alpha)
%amax: High bound for step
%x0: Current x value (initial x value if first iteration in an optimization
%script)
%c1: Constant parameter (0,1) typically between 0.5 and 0.95
%c2: Constant parameter (0,1) typically between 0.5 and 0.95

%Note: (c1,c2) = (0.5,0.5) leads to sufficient decrease in step size most
%of the time

%% OUTPUTS
%astar1: The "best" step size required to move to the next value of x, for
%the current iteration of an optimization algorithm

function [ astar1 ] = one_d_ls(amax,a1,x0,c1,c2)

%define objective function f here
f = @(x) -100*exp(-0.5*x)*sin(-2*x);
%gradient
g = @(x) -100*exp(-0.5*x)*(-0.5*sin(-2*x)-2*cos(-2*x));
%1D Hessian (2nd derivative)
B = @(x) 50*exp(-0.5*x)*((15/2)*sin(-2*x)-4*cos(-2*x));

%Evaluate at x_k
fk = f(x0); %equivalent to phi(0)
gk = g(x0); %equivalent to phi'(0)
Bk = B(x0);

%Find step direction
pk = -(1/Bk)*gk;

%Define scope of possible step sizes to search across
as = linspace(a1,amax,100);
as = [0 as];
i = 2;
ai = as(i);
while as(i) < amax
    %Define and evaluate univariate function phi and derivative phi'
    phi = f(x0+ai*pk);
    phiprime = g(x0+ai*pk)*pk;
    if (phi > (fk+c1*ai*gk*pk) || ( phi >= f(x0+as(i-1)*pk)) && i > 1)
        astar = extrapolate(as(i-1),ai);
        break
    end
    if abs(phiprime) <= (-c2*gk*pk)
        astar = ai;
        break
    end
    if phiprime >= 0 && i > 1
        astar = extrapolate(ai,as(i-1));
        break
    end
    i = i + 1;
    ai = as(i);
end
astar1 = astar;


%% Bisection Function to Interpolate Between Step Size Bounds
function [output1] = extrapolate(alo,ahi)
output = 0;
for k = 1:5
aj = (ahi-alo)/2;
phiaj = f(x0+aj*pk);
phialo = f(x0+alo*pk);
phiprime = g(x0+aj*pk)*pk;
if phiaj > (fk+c1*aj*gk) || phiaj >= phialo
    ahi = aj;
else
    if abs(phiprime) <= -c2*gk*pk
        output = aj;
        break
    end
    if phiprime*(ahi-alo) >= 0
        ahi = alo;
    end
    alo = aj;
end % if end
end % while end
output1 = output;
end %nested function end

end %function end


        