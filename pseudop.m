function [Mp] = pseudop(p)
% This function calculates gas viscosity (cP), compresivility factor (Z) 
% and pseudo-pressures (psi^2/cP) whit pressure data of the well (psi).
%
% INPUT DATA
% p             pressure data as a vector (psi)
% 
% OUTPUT DATA
% Mp            pseudo-pressures data as a vector (psi^2/cP) 
%
% Z & Visc_g will be calculated by a polynomial approximation whit pressure 
% data at GOR 62000 stf/STB.

% Z Coefficients
Z_a = -4.3114190846279e-16;
Z_b = 2.31875203787692e-12;
Z_c = 2.09683880523599e-8;
Z_d = -0.000109053802161887;
Z_e = 0.999688550578287;

% Z Polynomial
Z = Z_a*p.^4 + Z_b*p.^3 + Z_c*p.^2 + Z_d*p + Z_e;

% Gas Viscosity Coefficients
Visc_a = -1.02920472678361e-13;
Visc_b = 1.20448854132981e-9;
Visc_c = 1.52436143934205e-7;
Visc_d = 0.013447834921584;

% Gas Viscosity Polynomial
Visc_g = Visc_a*p.^3 + Visc_b*p.^2 + Visc_c*p + Visc_d;

% Check to see if pressure data is increasing or decreasing.
% m(p), Pseudo-pressures calculation using trapezoidal method.
if p(2) > p(1)
    Fun = 2*p./(Visc_g.*Z);
    mp = zeros(length(Fun),1);
    mp(1) = Fun(1)*p(1)/2;
    for i = 2:length(Fun)
        mp(i) = (Fun(i) + Fun(i-1))*(p(i)-p(i-1))/2;
    end
    Mp = cumsum(mp);
else
    p = flipud(p);
    Z = flipud(Z);
    Visc_g = flipud(Visc_g);
    Fun = 2*p./(Visc_g.*Z);
    mp = zeros(length(Fun),1);
    mp(1) = Fun(1)*p(1)/2;
    for i = 2:length(Fun)
        mp(i) = (Fun(i) + Fun(i-1))*(p(i)-p(i-1))/2;
    end
    Mp = flipud(cumsum(mp));
end
end
% REFERENCES
%
% https://www.sciencedirect.com/topics/engineering/pseudopressure.
% NARANJO AGUDELO, A. (2009). Evaluación de yacimientos de hidrocarburos. 
% Editorial Universidad Nacional de Colombia, sede Medellín.