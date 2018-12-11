%% Orbit Tools: Lambert Solver (Non-Universal, Elliptic Transfers Only)
%Benjamin Grace
function [orbittype,Vel,Vel_alt,Sol_info,t_p,min_en,min_ecc,V_ecc] = lambsolve(r1,r2,t,mu)
%% Inputs
%r1 = position vector of starting point 3x1 COLUMN VECTOR
%r2 = position vector of end point 3x1 COLUMN VECTOR
%t = transit time (t2-t1 in seconds)
%mu = gravitational parameter
%% Outputs
%orbittype = type of orbit: elliptic, hyperbolic, parabolic
%Vel = velocity vector set [V1 V2] for 0 <= theta < pi
%Vel_alt = velocity vector set [V1 V2] for pi <= theta < 2*pi
%Sol_info = further information about solutions including eccentricity and
%semimajor axis: [a; e; e_alt]
%t_p = parabolic transfer time
%min_en = minimum energy ellipse information: [a_m; e_m; t_m; t_m2]
%min_ecc = minimum eccentricity ellipse information: [a_ecc; e_ecc; t_ecc; t_ecc2]
%V_ecc = [V_ecc1 V_ecc2] velocity information for minimum eccentricity
%ellipse
%% Execution
%determine chord length: c
c = r2-r1;
%determine semiperimeter: s
s = 0.5*(norm(r1)+norm(r2)+norm(c));

%1) Compute limits on transfer:
%compute minimum semimajor axis required for an elliptic transfer (minimum
%energy ellipse)
a_m = (norm(r1)+norm(r2)+norm(c))/4;
%compute time on the minimum energy ellipse t_m
alpha_m = pi;
beta_m = 2*asin(sqrt(((s-norm(c))/s)));
t_m = (sqrt(s^3/8)*(alpha_m-beta_m+sin(beta_m)))/sqrt(mu);
alpha_m = 2*pi-alpha_m;
beta_m = -beta_m;
t_m2 = (sqrt(s^3/8)*(alpha_m-beta_m+sin(beta_m)))/sqrt(mu);
%compute orbital parameter for minimum energy ellipse
p_m = ((4*a_m*(s-norm(r1))*(s-norm(r2)))/norm(c)^2)*sin((alpha_m+beta_m)/2)^2;
e_m = sqrt(1-(p_m/a_m));
%compute minimum eccentricity for mininimum eccentricity ellipse
e_ecc = (norm(r2)-norm(r1))/norm(c);
%compute semimajor axis for minimum eccentricity ellipse
a_ecc = (norm(r1)+norm(r2))/2;

%2) Determine if transfer orbit is elliptical or hyperbolic:
%compute angle between vectors
theta = acos((r1'*r2)/(norm(r1)*norm(r2)));
%compute parabolic transfer time tp
if sin(theta) > 0
    sgn = 1;
end
if sin(theta) < 0
    sgn = -1;
end
t_p = ((sqrt(2)/3)*(s^(3/2)-sgn*(s-norm(c))^(3/2)))/sqrt(mu);
if t > t_p
    orbittype = "Elliptic";
else
    orbittype = "Hyperbolic";
end
if t == t_p
    orbittype = "Parabolic";
end
if strcmp(orbittype,"Elliptic") == 1

%3) Iteratively solve for semimajor axis a
    C = norm(c);
    %Seed initial guess
    ainit = a_m:0.01:500*a_m;
    for i = 1:length(ainit)
        finit(i) = (ainit(i)^(3/2))*(2*asin(sqrt(s/(2*ainit(i))))-2*asin(sqrt((s-C)/(2*ainit(i))))-sin(2*asin(sqrt(s/(2*ainit(i)))))+sin(2*asin(sqrt((s-C)/(2*ainit(i))))))-sqrt(mu)*t;
    end
    minval = min(abs(finit));
    initguess = ainit(find(abs(finit)-minval < 1E-10));
    a = initguess;
    %Initiate Newton Method
    while abs((a^(3/2))*(2*asin(sqrt(s/(2*a)))-2*asin(sqrt((s-C)/(2*a)))-sin(2*asin(sqrt(s/(2*a))))+sin(2*asin(sqrt((s-C)/(2*a)))))-sqrt(mu)*t) >= 1E-8
        func = (a^(3/2))*(2*asin(sqrt(s/(2*a)))-2*asin(sqrt((s-C)/(2*a)))-sin(2*asin(sqrt(s/(2*a))))+sin(2*asin(sqrt((s-C)/(2*a)))))-sqrt(mu)*t;
        C1 = -(C-s)*(2*a+C-s)/a^2;
        C2 = (2*a-s)*s/a^2;
        C3 = s/a;
        C4 = (-C+s)/a;
        dfunc = (a^2/(2*sqrt(a^3)))*(3*sqrt(C1)-3*sqrt(C2)+2*((C*(C-s)/sqrt(C1))-s^2/sqrt(C2)+s*(-C+s)/sqrt(C1))/a^2+6*asin(sqrt(C3)/sqrt(2))-6*asin(sqrt(C4)/sqrt(2)));
        a = a-(func/dfunc);
        if a < 0
            a = abs(a);
        end
        if a == a_m
            a = a_m;
            break
        end
    end

%4) Determine all possible solutions for alpha and beta:
a0 = a;
%for 0 <= theta < pi
beta01 = 2*asin(sqrt((s-norm(c))/(2*a0)));
%for pi <= theta < 2pi
beta02 = -beta01;
%shorter transfer time
alpha01 = 2*asin(sqrt(s/(2*a0)));
%longer transfer time
alpha02 = 2*pi-alpha01;
%compute parameters for terminal velocity equation
Ashort = sqrt(mu/(4*a0))*cot(alpha01/2);
Along = sqrt(mu/(4*a0))*cot(alpha02/2);
B1 = sqrt(mu/(4*a0))*cot(beta01/2);
B2 = sqrt(mu/(4*a0))*cot(beta02/2);

%5) Compute the terminal velocity vectors for both transfer angles,
%depending on transfer time given
%Note: sm corresponds to the small transfer angle, lg corresponds to the large transfer angle
if t > t_m
    V1_sm = (B1+Along)*(c./norm(c))+(B1-Along)*(r1./norm(r1));
    V2_sm = (B1+Along)*(c./norm(c))-(B1-Along)*(r2./norm(r2));
    V1_lg = (B2+Along)*(c./norm(c))+(B2-Along)*(r1./norm(r1));
    V2_lg = (B2+Along)*(c./norm(c))-(B2-Along)*(r2./norm(r2));
end
if t <= t_m
    V1_sm = (B1+Ashort)*(c./norm(c))+(B1-Ashort)*(r1./norm(r1));
    V2_sm = (B1+Ashort)*(c./norm(c))-(B1-Ashort)*(r2./norm(r2));
    V1_lg = (B2+Ashort)*(c./norm(c))+(B2-Ashort)*(r1./norm(r1));
    V2_lg = (B2+Ashort)*(c./norm(c))-(B2-Ashort)*(r2./norm(r2));
end

%6) Compute eccentricity information about each solution
if t > t_m
    p1 = ((4*a0*(s-norm(r1))*(s-norm(r2)))/norm(c)^2)*sin((alpha02+beta01)/2)^2;
    p2 = ((4*a0*(s-norm(r1))*(s-norm(r2)))/norm(c)^2)*sin((alpha02+beta02)/2)^2;
end
if t <= t_m
    p1 = ((4*a0*(s-norm(r1))*(s-norm(r2)))/norm(c)^2)*sin((alpha01+beta01)/2)^2;
    p2 = ((4*a0*(s-norm(r1))*(s-norm(r2)))/norm(c)^2)*sin((alpha01+beta02)/2)^2;
end
e1 = sqrt(1-(p1/a0));
e2 = sqrt(1-(p2/a0));

%7) Output information/solutions:
Vel = [V1_sm V2_sm];
Vel_alt = [V1_lg V2_lg];
Sol_info = [a0; e1; e2];
%minimum energy ellipse information: [a_m; e_m; t_m; t_m2]
min_en = [a_m; e_m; t_m; t_m2];
%minimum eccentricity ellipse information: [a_ecc; e_ecc; t_ecc; t_ecc2]
t_ecc = ((a_ecc^(3/2))*(2*asin(sqrt(s/(2*a_ecc)))-2*asin(sqrt((s-C)/(2*a_ecc)))-sin(2*asin(sqrt(s/(2*a_ecc))))+sin(2*asin(sqrt((s-C)/(2*a_ecc))))))/sqrt(mu);
t_ecc2 = ((a_ecc^(3/2))*(2*pi-2*asin(sqrt(s/(2*a_ecc)))+2*asin(sqrt((s-C)/(2*a_ecc)))-(sin(2*pi-2*asin(sqrt(s/(2*a_ecc)))))-sin(2*asin(sqrt((s-C)/(2*a_ecc))))))/sqrt(mu);
min_ecc = [a_ecc; e_ecc; t_ecc; t_ecc2];
%minimum eccentricity ellipse velocities
alphae = 2*asin(sqrt(s/(2*a_ecc)));
alphae2 = 2*pi-alphae;
betae = 2*asin(sqrt((s-norm(c))/(2*a_ecc)));
betae2 = -betae;
Ashorte = sqrt(mu/(4*a_ecc))*cot(alphae/2);
Alonge = sqrt(mu/(4*a_ecc))*cot(alphae2/2);
B1e = sqrt(mu/(4*a_ecc))*cot(betae/2);
B2e = sqrt(mu/(4*a_ecc))*cot(betae2/2);
V1_sm_ecc = (B1e+Ashorte)*(c./norm(c))+(B1e-Ashorte)*(r1./norm(r1));
V2_sm_ecc = (B1e+Ashorte)*(c./norm(c))-(B1e-Ashorte)*(r2./norm(r2));
V1_lg_ecc = (B2e+Alonge)*(c./norm(c))+(B2e-Alonge)*(r1./norm(r1));
V2_lg_ecc = (B2e+Alonge)*(c./norm(c))-(B2e-Alonge)*(r2./norm(r2));
V_ecc = [V1_sm_ecc V2_sm_ecc V1_lg_ecc V2_lg_ecc];
else
    Vel = [zeros(3,1) zeros(3,1)];
    Vel_alt = [zeros(3,1) zeros(3,1)];
    Sol_info = zeros(3,1);
    min_en = zeros(3,1);
    a_ecc = 0;
    t_ecc = 0;
    min_ecc = [a_ecc; e_ecc; t_ecc];
end
end