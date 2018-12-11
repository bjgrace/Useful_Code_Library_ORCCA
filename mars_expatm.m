%% Orbit Tools: Mars Exponential Atmospheric Model
%(accurate for 120-300 km)
%Benjamin Grace
function [density] = mars_expatm(h)
%% Inputs:
%r - altitude above planetary surface (km)

%% Outputs:
%density - density in kg/km^3

%% Execution:
%check for altititude range
if h >= 120 && h < 200
    C1 = -2.55314E-10;
    C2 = 2.31927E-07;
    C3 = -8.33206E-05;
    C4 = 0.0151947;
    C5 = -1.52799;
    C6 = 48.69659;
    density = (1E9)*exp(C1*(h^5)+C2*(h^4)+C3*(h^3)+C4*(h^2)+C5*(h)+C6);
end
if h >= 200 && h <= 300
    C1 = 2.65472E-11;
    C2 = -2.45558E-08;
    C3 = 6.31410E-06;
    C4 = 4.73359E-04;
    C5 = -0.443712;
    C6 = 23.79408;
    density = (1E9)*exp(C1*(h^5)+C2*(h^4)+C3*(h^3)+C4*(h^2)+C5*(h)+C6);
end

end