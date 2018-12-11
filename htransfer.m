%% Orbit Tools: Hohmann Transfer (Type I, Type II, etc.)
%NOT FOR INTERPLANETARY USE (yet)
function [deltaV1, deltaV2] = htransfer(r0,a0,rdes,ades,mu,type)
%% Inputs:
%r0 - current radius (m or km)
%a0 - current semi-major axis (m or km) (same as r if on circular orbit)
%rdes - desired radius (rp or ra) (m or km)
%ades - desired semimajor axis (if target orbit is not circular), otherwise
%can just plug in r, or an arbitrary value and the script will ignore it)
%(m or km)
%mu - gravitational parameter (units corresponding to r,a)
%type - "circular" or "elliptic" determines the type of the target orbit

%% Outputs:
%deltaV1 - the required velocity change for the first burn
%deltaV2 - the required velocity change for the second burn

%% Execution:
%magnitude of all quantities:
r0 = norm(r0);
a0 = norm(a0);
rdes = norm(rdes);
ades = norm(ades);
%determine current orbital velocity and transfer orbit semi-major axis
v0 = sqrt(2*(mu/r0)-mu/a0);
at = (r0+rdes)/2;

%for circular orbit
if strcmp(type,"circular") == 1
%speed at desired circular orbit
vdes = sqrt(mu/rdes);
%speed of transfer orbit (at current radius)
vt1 = sqrt(2*(mu/r0)-mu/at);
%speed of transfer orbit (at desired radius)
vt2 = sqrt(2*(mu/rdes)-mu/at);
deltaV1 = vt1-v0;
deltaV2 = vdes-vt2;
end

%for elliptic orbit
if strcmp(type,"elliptic") == 1
%speed at desired elliptic orbit
vdes = sqrt(2*(mu/rdes)-mu/ades);
%speed of transfer orbit (at current radius)
vt1 = sqrt(2*(mu/r0)-mu/at);
%speed of transfer orbit (at desired radius)
vt2 = sqrt(2*(mu/rdes)-mu/at);
deltaV1 = vt1-v0;
deltaV2 = vt2-vt1;
end
end
