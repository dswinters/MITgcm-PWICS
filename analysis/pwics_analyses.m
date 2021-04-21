%% Example script for running batch analyses

clear all, close all

% Add some paths for external code (may need to be modified)
addpath ../setup % stratification code used in gendata
addpath ~/MATLAB/gsw % contains inertial freq code
addpath ~/MITgcm/MITgcm/utils/matlab % MITgcm data handling code
addpath ~/MATLAB % colormaps - I used 'cmocean', which I keep here

% theta and kTopo setup copied from gendata.m
% TODO: Reorganize input generation so that these are only defined in one place
theta = [0 20 40 60]; % angle of PW wave vector wrt line perpendicular to coast [deg]
lTopo = linspace(20e3,700e3,30);
kTopo = 2*pi./lTopo;

for i = 1:length(theta)
    for j = 1:length(kTopo)
        calcFluxBC(theta(i),kTopo(j));
    end
end

plotReson(theta,kTopo)

load plotReson.mat;
for i = 1:length(theta)
    [~,idx] = max(aflx{i});
    plotVel(theta(i),kTopo(idx))
end

$./plotVel_movie.sh
