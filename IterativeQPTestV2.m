%% Iterative Quadratic Program - Inequality Constraints
%% Active Set QP Algorithm
%% Source:
%Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (Second ed.). Springer. 
%doi:10.1007/978-0-387-40065-5
%Chapter 16: Quadratic Programming - Active Set Method for Convex QP's

%% Instructions for use:
%{
-FOR USE WITH WELL DEFINED QP WITH LINEAR CONSTRAINTS ONLY
-Better off using the Active Set SQP to solve both QP and SQP with active
set method, this script is fairly rudamentary and requires more manual
input
-Alter appropriate variables as listed below
-Make a feasible guess for a working set
(Aw,Bw such that A*x-b <= 0, if A*x-b > 0, constraint is not in the working set)
-For problem formulation min(x) q(x) = x'*G*x + c, s.t. A*x >= b
%}

%% INPUTS
%{
G: Hessian, G, from q(x) = x'*G*x + c
c: Gradient, c, from q(x) = x'*G*x + c
x0: Initial state vector
b: Constraint vector in A*x >= b
A: Constraint matrix in A*x >= b
Aw: Feasible working set, constraint matrix (i.e. ([A(1,:); A(4,:)])
bw: Feasible working set, constraint vector (i.e. ([b(1,:); b(4,:)])
%}
%% OUTPUTS
%{
xfinal: Final value of states at minimum
%}

tic
%coeff = [1 0 0 1 -2 -5]; %coefficients for quadratic problem a11x1x1, a12x1x2, ... etc.
const = [1 4.25]; %constants in problem, independent of x
G = [2 0; 0 2];%[4 2; 0 1];%[1 -2; 0 2];%[4 2; 0 1];%[2 0; 0 2];
c = [-2 -5];%[2 3];%[-2 -6];%[2 3];%[-2 -5];
x0 = [0; 0;]; %initial x vector
v = length(x0); %number of state variables
%coeff(length(coeff)-v+1:length(coeff)); %costs (constants in front of linear terms)
b = [-2; -6; -2; 0; 0];%[0; -4; -3];%[-1; -2; 0; 0];%[0; -4; -3];%[-2; -6; -2; 0; 0]; %constraints
A = [1 -2; -1 -2; -1 2; 1 0; 0 1];%[1 -1; -1 -1; -1 0];%[-0.5 -0.5; 1 -2; 1 0; 0 1];%[1 -1; -1 -1; -1 0];%[1 -2; -1 -2; -1 2; 1 0; 0 1]; %coefficients for constraints
k = 1;
%find working set
Aw = [A(1,:); A(4,:);];
bw = [b(1,:); b(4,:);];
W0 = 1:1:length(Aw(:,1));
exit = 0;
while k > 0
    %{
    if k == 1
        G = HessianApprox(v,x0,coeff,const);
    else
        G = G;
    end
    %}
    g = c'+G*x0;
    h = Aw*x0-bw;
    K = [G Aw'; Aw zeros(length(Aw(:,1)))];
    plam = K\[g; h];
    xn = x0 + -plam(1:v);
    if isnan(xn) > 0
        disp('Infeasible Working Set!');
        break
    end
    %yk = JacobianGradApprox(v,xn,coeff,const)'-JacobianGradApprox(v,x0,coeff,const)';
    %sk = xn-x0;
    %B = G;
    k = k+1;
    xstore(:,k) = xn;
    %G = B+(((yk-B*sk)*(yk-B*sk)')/((yk-B*sk)'*sk));
    neg = find(plam(v+1:length(plam))<0);
    size = numel(neg);
    repeat = 0;
    if size > 0
        if length(Aw(:,1)) == 1
            Aw(:,:) = [];
            bw(:,:) = [];
        else
        Aw = Aw(1:length(W0)-1,:);
        bw = bw(1:length(W0)-1,:);
        end
    else
    pn = xn-x0;
    count = 1;
    while count <= length(A(:,1))
        if A(count,:)*xn-b(count,:) < 0
            xt = A(count,:)\b(count,:);
            alpha = pn\(xt-xn);
            xn = xn+alpha*pn;
            Aw = A(count,:);
            bw = b(count,:);
            xstore2(:,k) = xn;
            repeat = 1;
        else
            exit = 1;
        end
        count = count+1;
    end
    if exit == 1 && repeat ~= 1
        break
    end
    end
end

    % compute p
    %if 0, check lagrange multipliers
    % if all positive, return x
    %otherwise remove one negative lagrange multiplier from wk+1
    %update xk+1
    %end if
    %computer ak from xk+1 = xk+akpk
    %if encounter blocking constraint add to wk+1
    %else
    %wk+1 = wk
toc

xfinal = xn; %display the new/final state

%% Plotting Script (uncomment to use)
%{
% generate function space
[X1,X2] = meshgrid(0:0.1:2,0:0.1:2);
Zf = (X1-ones(length(X1(:,1)),length(X1(1,:)))).^2+(X2-2.5*ones(length(X2(:,1)),length(X2(1,:)))).^2;

X1C = [1:0.1:4];
X2C1 = (-X1C-2*ones(length(X1C),1))./-2;
X2C2 = (X1C-6*ones(length(X1C),1))./-2;
X2C3 = (X1C-2*ones(length(X1C),1))./2;

%compute constraints
fc1 = @(x1) 0.5*x1+1;
fc2 = @(x1) -0.5*x1+3;
fc3 = @(x1) 0.5*x1-1;
fc4 = @(x1) 0;
fc5 = @(x1) 0;
f1 = @(x1) 1;

%evaluate and plot contours of function space
figure(1)
mesh(X1,X2,Zf); hold on
%fplot(X1C,X2C1,X1C,X2C2,X1C,X2C3); hold on
fplot(fc1); hold on;
fplot(fc2); hold on;
fplot(fc3); hold on;
fplot(fc4); hold on;
fplot(fc5); hold on;
plot3([0 0 0 0 0],[0 1 2 3 4],[0 2 4 6 8]); hold on
plot3([0 1 2 3 4],[0 0 0 0 0],[0 2 4 6 8]); hold on
plot3(1.4,1.7,0.8,'-o');
colorbar
title('Constrained Minimum of f(x) - Active Set Method QP');
xlabel('X1');
ylabel('X2');
zlabel('f(X1,X2)');
xlim([-2 4]);
ylim([-1 3]);

%evaluate and plot iterations of optimization 
figure(2)
plot(xstore(1,:),xstore(2,:),'-o'); hold on
fplot(fc1); hold on;
fplot(fc2); hold on;
fplot(fc3); hold on;
fplot(fc4); hold on;
fplot(fc5); hold on;
fplot(f1); hold on;
xlim([0 5])
ylim([0 5])
title('Iterations of X* Constraints 3,5');
xlabel('X1');
ylabel('X2');
grid on
%}