% Use this function to get setup parameters to avoid hard-coding parameters in
% multiple places.

function params = gendata_params()

params.theta = [-60:20:60]; % angle of PW wave vector wrt line perpendicular to coast [deg]
params.theta = [20]; % angle of PW wave vector wrt line perpendicular to coast [deg]
params.lprof = 20e3; % depth (negative) or extent (positive) of corrugations
params.lTopo = linspace(20e3,700e3,30);
params.kTopo = 2*pi./lTopo;
params.f = 1e-4;
params.g = 9.81;
params.om = 1.36*params.f;
params.deltaT = 500;
