%% Orbit Tools: Numerical Integration Package
%Benjamin Grace

function [time,r,rdot,elements,energy,momentum] = orbintegrator(r0,rdot0,mu,ti,tf,dt,u,reltol,abstol)
%% Inputs
%r0 - initial position vector (in m or km, 3x1 column vector)
%rdot0 - initial velocity vector (in m/s or km/s, 3x1 column vector)
%mu - gravitational parameter in units corresponding to r, rdot
%ti - initial time
%tf - final time
%dt - step size
%u - perturbation (3x1 column vector)
%reltol - relative tolerance for ode45 1e-12 recommended
%abstol - absolute tolerance for ode45 1e-14 recommended
%% Outputs
%time - an Nx1 time vector with all timesteps from 0 to tf
%r - a 3xN matrix of position vectors corresponding to time
%rdot - a 3xN matrix of velocity vectors corresponding to time
%elements - a 6xN matrix of orbital elements vectors corresponding to time [a,e,i,RAAN,w,f]
%energy - an Nx1 vector of orbital energy
%momentum- a 3xN vector of orbital angular momentum

%% Execution
%numerically integrate orbital equation for 2 body problem
%tvec = 0:dt:tf;
opts = odeset('RelTol',reltol,'AbsTol',abstol);
%r = zeros(3,length(tvec));
%rdot = zeros(3,length(tvec));
%r(:,1) = r0;
%rdot(:,1) = rdot0;
%for i = 2:length(tvec)
tvec = [ti:dt:tf];
if dt ~= 0
    [syst,sysdata] = ode45(@(t,y) twobodyODE(t,y,mu,u),tvec,[r0; rdot0],opts);
end
if dt == 0
    [syst,sysdata] = ode45(@(t,y) twobodyODE(t,y,mu,u),[ti tf],[r0; rdot0],opts);
end
    r = sysdata(:,1:3)';
    rdot = sysdata(:,4:6)';
%end
time = syst;
energy = zeros(length(time),1);
momentum = zeros(3,length(time));
%compute orbit energy and momentum at each timestep
for i = 1:length(time)
    energy(i) = (1/2)*norm(rdot(:,i))^2-mu/norm(r(:,i));
    momentum(:,i) = cross(r(:,i),rdot(:,i));
end
%evaluate orbit elements at each timestep
%preallocate
elements = zeros(6,length(time));
for i = 1:length(time)
    evec = cross((1/mu)*rdot(:,i),momentum(:,i))-(r(:,i)/norm(r(:,i)));
    %compute eccentricity and momentum scalars
    e = norm(evec);
    h = norm(momentum(:,i));
    %compute semimajor axis from orbit parameter
    p = h^2/mu;
    a = p/(1-e^2);
    %compute inclination
    hnorm = momentum(:,i)./h;
    inc = acos(dot(hnorm,[0;0;1]));%norm(atan2(cross(hnorm,[0; 0; 1]),[0 0 1]*hnorm));
    %compute node vector
    node = cross([0;0;1],hnorm)./norm(cross([0;0;1],hnorm));
    RAAN = atan2([0 1 0]*node,[1 0 0]*node);
    if isnan(RAAN) == 0 && RAAN < 0
        RAAN = RAAN+2*pi;
    end
    %compute perpendicular node vector
    nodeperp = cross(hnorm,node);
    %compute argument of periapsis
    w = atan2(evec'*nodeperp,evec'*node);
    if isnan(w) == 0 && w < 0
        w = w+2*pi;
    end
    %compute true anomaly
    eperp = cross(hnorm,evec./e);
    f = atan2(r(:,i)'*eperp,r(:,i)'*(evec./e));
    elements(1,i) = a;
    elements(2,i) = e;
    elements(3,i) = rad2deg(inc);
    elements(4,i) = rad2deg(RAAN);
    elements(5,i) = rad2deg(w);
    elements(6,i) = rad2deg(f);
    if isnan(elements(4,i)) == 1
        elements(4,i) = 0;
    end
    if isnan(elements(5,i)) == 1
        elements(5,i) = 0;
    end
    if isnan(elements(6,i)) == 1
        elements(6,i) = 0;
    end
end

    function [dydt] = twobodyODE(t,y,mu,u)
        dydt = zeros(6,1);
        dydt(1) = y(4);
        dydt(2) = y(5);
        dydt(3) = y(6);
        dydt(4) = -mu*y(1)/(sqrt(y(1)^2+y(2)^2+y(3)^2))^3+u(1);
        dydt(5) = -mu*y(2)/(sqrt(y(1)^2+y(2)^2+y(3)^2))^3+u(2);
        dydt(6) = -mu*y(3)/(sqrt(y(1)^2+y(2)^2+y(3)^2))^3+u(3);
    end
end