%% PQP - Parallel Quadratic Program
%% Notes:
%{
-For use with linear MPC
-This is framework for the PQP used in the Simulink Linear MPC model
-Works as long as inputs are specified, is not in function form like the
PQP in the Simulink model
%}
%% Source
%{
Brand, M., Shilpiekandula, V., Yao, C., & Bortoff, S. A. (2011). A Parallel Quadratic 
Programming Algorithm for Model Predictive Control. Mitsubishi Electric 
Research Laboratories, TR2011(056). Retrieved from http://www.merl.com
%}

%% INPUTS/ Changeable Variables
%{
itermax: max number of iterations (define in script)
e: tolerance (define in script)
Ad: Discretized state transfer matrix (A) from DLQR (initialize beforehand)
Bd: Discretized input transfer matrix (B) from DLQR (initialize beforehand)
Q: Quadratic weights from objective function (initialize beforehand)
R: Resource weights from objective function (initialize beforehand)
S: Terminal weights from objective function (initialize beforehand)
x0: Initial state (initialize beforehand)
xN: Final/Target state (initialize beforehand)
can define u0 as initial input, not needed since sizes gained from Ad and Bd
umin: 1*length(u0) vector of minimum input bounds (define in script)
umax: 1*length(u0) vector of maximum input bounds (define in script)
xmin: 1*length(x0) vector of minimum state bounds (define in script)
xmax: 1*length(x0) vector of maximum state bounds (define in script)
%}

%% OUTPUTS
%{
Ustar: Calculated control input
%}

% Initialization (matrix dimensions, steps, state space)
itermax = 10000;
e = 10; %tolerance
n = length(Ad(:,1));
u = length(Bd(1,:));
Qlarge = zeros(N*u+(N)*n);
QN = zeros(N*u+(N)*n);
i = 1;
for ncount = 1:N
    QN(i:(ncount*u),i:(ncount*u)) = R;
    i = i+u;
end
for ncount = 1:(N-1)
    QN(i:(i+n-1),i:(i+n-1)) = Q;
    i = i+n;
end
QN(i:(i+(n-1)),i:(i+(n-1))) = S;
    i = i+n;
    
psi = zeros(n*N,n);
omega = zeros(n*N,u*N);
L1 = zeros(u*N,u*N);
L2 = zeros(n*N,n*N);
L4 = zeros(n,n*N);
L1 = QN(1:u*N,1:u*N);
L2 = QN(u*N+1:(u*N+n*N),u*N+1:(u*N+n*N));
L4(1:n,length(L4(1,:))-n+1:length(L4(1,:))) = eye(n);

storem = zeros(n,n,N);
i = 1;
while i <= N
    storem(:,:,i) = Ad^(i);
    i = i+1;
end
ncount = 1;
i = 1;
for ncount = 1:N
    psi(i:ncount*n,1:n) = storem(:,:,ncount);
    i = i+n;
end

storemat2 = zeros(n,u,N);
i = 1;
while i <= N
    storemat2(:,:,i) = Ad^(i-1)*Bd;
    i = i+1;
end
ncount = 1;
N2 = 1;
j = 1;
i = 1;
skipd = 0;
subd = 0;
while N2 <= N
while ncount <= N-subd
        omega((i+skipd):(skipd+ncount*n),j:N2*u) = storemat2(:,:,(ncount));
        i = i+n;
        ncount = ncount+1;
end
i = 1;
ncount = 1;
j = j+u;
N2 = N2+1;
skipd = skipd+n;
subd = subd+1;
end

x0 = [0; 0.5; 0; 2; 2; 2];
xN = [0.7; 0.7; 0.5; 3; 3; 3];
iniy = ones(n*N,1);

% Set minimums and maximums
Xkmax = zeros(1,n*N);
Xkmin = zeros(1,n*N);
Ukmax = zeros(1,u*N);
Ukmin = zeros(1,u*N);
umax = [0.2 0.2 0.2];
umin = [-0.2 -0.2 -0.2];
xmax = [1 1 1 2 2 2];
xmin = [-1 -1 -1 -2 -2 -2];
i = 1;
while i <= u*N-u+1
    Ukmax(i:i+u-1) = umax;
    Ukmin(i:i+u-1) = umin;
    i = i+1;
end
i = 1;
while i <= n*N-n+1
    Xkmax(i:i+n-1) = xmax; 
    Xkmin(i:i+n-1) = xmin;
    i = i+1;
end

%Form matrices necessary for primal and dual functions

Qp = 2*L1+2*omega'*L2*omega;
Hp = (2*x0'*psi'*L2*omega-2*xN'*S*L4*omega)';
V = [eye(u*N); -eye(u*N)];
W = [Ukmax, -Ukmin]';
y = iniy;
ydiag = zeros(n*N,n*N);
%initialize diag(y);
i = 1;
j = 1;
while i <= n*N
    while j <= n*N
        if i == j
            ydiag(i,j) = y(i);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end

%Solve dual problem
Qinv = inv(Qp);
Qdual = V*Qinv*V';
hdual = W+V*Qinv*Hp;

Qplus = zeros(n*N,n*N);
Qminus = zeros(n*N,n*N);
hplus = zeros(n*N,1);
hminus = zeros(n*N,1);

%select suitable r value
Qminusmax = Qplus;
hplusdiag = Qplus;
i = 1;
j = 1;
while i <= n*N
    while j <= n*N
        Qminusmax(i,j) = (1/2)*(-Qdual(i,j)+abs(-Qdual(i,j)));
        if i == j
            hplusdiag(i,j) = (1/2)*(hdual(i)+abs(hdual(i)));
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end

Knn_r = hplusdiag*inv(ydiag)+Qminusmax; %0-r == Knn-r
rvec = -max(Knn_r);
r = 0;%max(rvec);
rdiag = zeros(n*N);

i = 1;
j = 1;
while i <= n*N
    while j <= n*N
        if i == j
            rdiag(i,j) = r;
        end
    j = j+1;
    end
    j= 1;
    i = i+1;
end

i = 1;
j = 1;
while i <= n*N
    while j <= n*N
        Qplus(i,j) = (1/2)*(Qdual(i,j) + abs(Qdual(i,j))) + rdiag(i,j);
        Qminus(i,j) = (1/2)*(-Qdual(i,j) + abs(-Qdual(i,j))) + rdiag(i,j);
        if i == j
            hplus(i) = (1/2)*(hdual(i)+abs(hdual(i)));
            hminus(i) = (1/2)*(-hdual(i)+abs(-hdual(i)));
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end

count = 1;
i = 1;
while count <= itermax
    while i <= n*N
        Qminy = Qminus*y;
        Qplusy = Qplus*y;
        y(i) = y(i)*((hminus(i)+Qminy(i))/(hplus(i)+Qplusy(i)));
        i = i+1;
    end
    Qminy = Qminus*y;
    Qplusy = Qplus*y;
    x = diag(y)*inv(diag(Qplusy+hplus))*(Qminy+hminus);
    if abs(x-y) <= ones(n*N,1)*e
        ystar = y;
        break
    end
    if count == itermax
        ystar = y;
    end
    i = 1;
    count = count+1;
end

%Solve primal problem
Ustar = -Qinv*(Hp+V'*ystar);


    