%% Orbit Tools: Orbital Elements to State Vector
%Benjamin Grace
function [rvec,rdotvec] = orbel2state(mu,a,e,i,RAAN,w,ft,opt)
%% Inputs
%NOTE: ALL ANGULAR QUANTITIES MUST BE IN RADIANS
%mu = gravitational parameter of central body
%a = semimajor axis (km)
%e = eccentricity
%i = inclination (rad)
%RAAN = right ascension (rad)
%w = argument of perigee (rad)
%ft = true anomaly (rad) OR time in orbit (if opts == "time") 
%opt = "time" or "anomaly" depending on known information

%% Outputs
%NOTE: Outputs will be in km or m depending on form of mu used
%rvec = position state vector in cartesian coordinates (km/s, m/s)
%rdotvec = velocity state vector in cartesian coordinates (km/s, m/s)

%% Execution
if strcmp(opt,"anomaly")==1
    f = ft;
    %Elliptical/Parabolic Case
    if e <= 1
        %compute orbital parameter p:
        p = a*(1-e^2);
        %compute scalar orbital momentum:
        h = sqrt(p*mu);
        %compute scalar r from conic equation:
        r = p/(1+e*cos(f));
        %compute scalar rdot from conic equation and momentum:
        rdot = (e*sin(f)*h)/p;
        %compute time derivative of true anomaly
        fdot = h/(r^2);
    end
    %Hyperbolic Case
    if e > 1
        %compute orbital parameter p:
        p = a*(e^2-1);
        %compute scalar orbital momentum:
        h = sqrt(p*mu);
        %compute scalar r from conic equation:
        r = p/(1+e*cos(f));
        %compute scalar rdot from conic equation and momentum:
        rdot = (e*sin(f)*h)/p;
        %compute time derivative of true anomaly
        fdot = h/(r^2);
    end   
end
if strcmp(opt,"time") == 1
    t = ft;
    %Elliptical/Parabolic Case
    if e <= 1
        %Solve Kepler's equation using time
        M = sqrt(mu/a^3)*t;
        E = M;
        funcE = M-E+e*sin(E);
        %begin newton solver
        while abs(funcE) >= 1E-8
            funcE = M-E+e*sin(E);
            dE = -1+e*cos(E);
            E = E-(funcE/dE);
        end
        %end newton solver
        f = 2*atan2(tan(E/2)*sqrt(1+e),sqrt(1-e));
        if f < 0
            f = f+2*pi;
        end
        %compute orbital parameter p:
        p = a*(1-e^2);
        %compute scalar orbital momentum:
        h = sqrt(p*mu);
        %compute scalar r from conic equation:
        r = p/(1+e*cos(f));
        %compute scalar rdot from conic equation and momentum:
        rdot = (e*sin(f)*h)/p;
        %compute time derivative of true anomaly
        fdot = h/(r^2);
    end
    %Hyperbolic Case
    if e > 1
        %Solve hyperbolic form of Kepler's equation using time
        M = sqrt(mu/abs(a)^3)*t;
        F = M;
        %begin newton solver
        funcF = M-F+e*sinh(F);
        while abs(funcF) >= 1E-8
            funcF = M-F+e*sinh(F);
            dF = -1+e*cosh(F);
            F = F-(funcF/dF);
        end
        %end newton solver
        f = 2*atan2(tanh(F/2)*sqrt(e+1),sqrt(e-1));
        %compute orbital parameter p:
        p = a*(e^2-1);
        %compute scalar orbital momentum:
        h = sqrt(p*mu);
        %compute scalar r from conic equation:
        r = p/(1+e*cos(f));
        %compute scalar rdot from conic equation and momentum:
        rdot = (e*sin(f)*h)/p;
        %compute time derivative of true anomaly
        fdot = h/(r^2);
    end
end
%compute r,v -> ECI coordinates transformation
l1 = cos(RAAN)*cos(w)-sin(RAAN)*sin(w)*cos(i);
l2 = -cos(RAAN)*sin(w)-sin(RAAN)*cos(w)*cos(i);
m1 = sin(RAAN)*cos(w)+cos(RAAN)*sin(w)*cos(i);
m2 = -sin(RAAN)*sin(w)+cos(RAAN)*cos(w)*cos(i);
n1 = sin(w)*sin(i);
n2 = cos(w)*sin(i);
zeta = r*cos(f); 
eta = r*sin(f);
zetadot = rdot*cos(f)-r*sin(f)*fdot;
etadot = rdot*sin(f)+r*cos(f)*fdot;
%compute ECI position and velocity
rvec = [l1*zeta+l2*eta; m1*zeta+m2*eta; n1*zeta+n2*eta];
rdotvec = [l1*zetadot+l2*etadot; m1*zetadot+m2*etadot; n1*zetadot+n2*etadot];




