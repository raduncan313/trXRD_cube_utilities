clear all
close all

% This is an example script for how to use the Cube class (defined in
% Cube.m) for analyzing time-resolved x-ray scattering experimental data.
% Please ask me (Ryan) for the .h5 and .vasp files needed to run the
% examples below.

%% Load a delay scan cube and inspect using the xdscode plotdatacube function
addpath(genpath('xdscode'));
cubedir = 'C:/Users/radun/OneDrive - Stanford/Documents/SLAC/Reis Group/TSI/xpply5120/cubes/';
t0 = 9; % Nominal time-zero for this dataset.
runs = {253}; % Run numbers in cell array (multiple runs are summed together into one cube)
scan_var_name = 'delay'; % Type of scan ('delay', 'phi', 'theta', 'chi', 'mu', etc.)
info = 'testing'; % info string for run -- used in plot titles/legends and saved filenames

d = load_cube(cubedir, runs, scan_var_name, info, t0); % Using helper function defined at end of script
f = d.plotdatacube('on', [-6 -1], 'log') % uses `plotdatacube` from `xdscode` library

%% Helper functions

function d = load_cube(cubedir, runs, scan_var_name, info, t0)
    d = Cube(cubedir, runs, scan_var_name, info);
    d.transpose();
    d.norm_i0();
    d.subtract_t0(t0);
end