%% Orbit Tools: State Vector to Orbital Elements
%Benjamin Grace
function [elements,t_periapsis] = state2orbel(r,rdot,mu,t0)
%% Inputs:
%r = position state in km
%rdot = velocity state in km/sec
%mu = gravitational parameter for central body
%t0 = current time (optional, input only if calculating time of periapsis
%passage)

%% Outputs:
%elements = 6x1 vector of orbital elements: [a e i RAAN w f]'
%^ the output is in degrees
%t_periapsis = time of periapsis passage (optional)

%% Execution
%allocate space for elements
elements = zeros(6,1);
%compute momentum and eccentricity vectors from state
hvec = cross(r,rdot);
evec = cross((1/mu)*rdot,hvec)-(r/norm(r));
%compute eccentricity and momentum scalars
e = norm(evec);
h = norm(hvec);
%compute semimajor axis from orbit parameter
p = h^2/mu;
a = p/(1-e^2);
%compute inclination
hnorm = hvec./h;
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
f = atan2(r'*eperp,r'*(evec./e));
elements(1) = a;
elements(2) = e;
elements(3) = rad2deg(inc);
elements(4) = rad2deg(RAAN);
elements(5) = rad2deg(w);
elements(6) = rad2deg(f);

if isnan(elements(4)) == 1
    elements(4) = 0;
end
if isnan(elements(5)) == 1
    elements(5) = 0;
end
if isnan(elements(6)) == 1
    elements(6) = 0;
end
%compute time of periapsis (if time is given)
switch nargin
    case 4
        if e < 1
            E = 2*atan2(sqrt(1-e)*tan(f/2),sqrt(1+e));
            t_periapsis = t0-sqrt(a^3/mu)*(E-e*sin(E));
        end
        if e > 1
            F = 2*atanh(sqrt((e-1)/(e+1))*tan(f/2));
            t_periapsis = t0+sqrt(abs(a)^3/mu)*(F-e*sinh(F));
        end
        if e == 1
            t_periapsis = t0-((1/2)*sqrt(p^3/mu)*tan(f/2)*(1+(1/3)*tan(f/2)^2));
        end
    case 3
        t_periapsis = 0;
end
end