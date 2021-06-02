%% Example script for running batch analyses

clear all, close all

% Add some paths for external code (may need to be modified)
addpath ../setup % stratification code used in gendata
addpath ~/MATLAB/gsw % contains inertial freq code
addpath ~/MITgcm/MITgcm/utils/matlab % MITgcm data handling code
addpath ~/MATLAB % colormaps - I used 'cmocean', which I keep here

cd ../setup
params = gendata_params();
theta = params.theta;
lTopo = params.lTopo;
kTopo = params.kTopo;
cd ../analysis

for i = 1:length(theta)
    for j = 1:length(kTopo)
        calcFluxBC(theta(i),kTopo(j));
    end
end

plotReson(theta,kTopo);

load plotReson.mat;
for i = 1:length(theta)
    [~,idx] = max(aflx{i});
    plotVel(theta(i),kTopo(idx))
    compute_ctw_wavelength(theta(i),kTopo(idx));
end

% !./plotVel_movie.sh
