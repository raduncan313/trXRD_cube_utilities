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

%% Apply a threshold mask, plot the total masked signal, and save it as a .csv file
th = 0.05;
f1 = d.thresh_and_plot(th); % applies a mask that zeros everything dimmer than 100*th % of the max integrated pixel value, and plots the resulting signal integrated over the whole detector
f2 = d.plot_mask(1); % show the mask applied in the previous line
folder = 'sigs/';
d.write_sigs_to_csv(folder); % Saves the signals from the mask as .csv files in `folder`.

%% Draw ROIs and plot the integrated signals
n = 1;
f3 = d.plot_rois(n) % draw n ROIs and plot the integrated signals

%% Apply LPSVD autoregressive modeling to the signal generated in the previous section
L = 8 % number of singular values to keep in the SVD decomposition
rat = 0.75; % ratio of the number of autoregression coefficients relative to the number of time steps
t_start = 3 % time point after which to begin the modeling
f4 = d.LPSVD_sig(1, L, rat, t_start)

%% Plot a dispersion relation by drawing a lineout on the detector image
f5 = d.plot_lineout();

%% Auto-detect ROIs using PCA-OCSVM-DBSCAN sequence
numcomponents = 10;
epsilon = 8;
minpts = 10;
frac = 0.03;
f6 = d.auto_signal_figure(numcomponents, epsilon, minpts, frac);

%% Auto-detect ROI with GUI

params = struct();
params.nc = 10;
params.frac = 0.03;
params.minpts = 10;
params.eps = 8;

f = d.auto_signal_GUI(params);

%% Construct an xdscode `geometry` struct, load an angle scan, and map it into reciprocal space

% The following several lines construct a geometry structure from the
% xdscode library, which will be used to configure a cube and perform hkl
% mapping

% read crystal structure info
[P, tau_p, attypes_p, atnums_p] = readPOSCAR('TaSe42I_sm_isp_sd_1502042_springermaterials.vasp');
tau_p = P*tau_p; % convert to cartesian
at = P;
bg = inv(at);
geometry.realvecs = at';
geometry.primvects = bg';
geometry.basis=tau_p';
geometry.lambda0 = 12.398/9.55; % x-ray wavelength

k0 = 1/geometry.lambda0*[1,0,0];

% Configure detector info
geometry.detector.det_dist = 600;
geometry.detector.det_pixels_horz=1030;
geometry.detector.det_pixels_vert=1064;
geometry.detector.det_size_horz = 0.075*geometry.detector.det_pixels_horz;
geometry.detector.det_size_vert = 0.075*geometry.detector.det_pixels_vert;
geometry.imageNy=round(geometry.detector.det_pixels_horz/10);
geometry.imageNz=round(geometry.detector.det_pixels_vert/10);
geometry.beam_center = [580, 534];

% Goniometer angles
geometry.theta = -0.5;
geometry.phi = -22.2;
geometry.chi = 90;
geometry.mu = 0;

% Detector angles
geometry.nu = -24.57;
geometry.delta = 20.6;

% Sample orientation
sample_normal=[1,1,0];
sample_inplane=[0,0,1];
geometry.SamRot = rotationmat3D(45,[0,0,1]);
geometry.rot_matrix = @(a,b,c,d) huber_sixcircle_matrix_lw89(a,b,c,d);

% Configure ROI to be mapped
ROI_lim_h = [300 550];
ROI_lim_v = [100 500];
ROI_lims = {ROI_lim_h, ROI_lim_v};

% pixel intensity lower and upper thresholds
threshs = [40 1000];

% Now configure the cube with the geometry structure included:

cubedir = './';
runs = {129};
scan_var_name = 'phi';
info = 'phiscan';
t0 = 0;

d2 = load_cube_phiscan(cubedir, runs, scan_var_name, info, geometry);
f6 = d2.hklmap(ROI_lims, threshs);

%% Helper functions

function d = load_cube(cubedir, runs, scan_var_name, info, t0)
    d = Cube(cubedir, runs, scan_var_name, info); % instantiate and construct Cube object
    d.transpose(); % transpose the detector dimensions if needed
    d.norm_i0(); % divide by i0 values (accounting for XFEL intensity fluctuation)
    d.subtract_t0(t0); % subtract t0 from the scan_var values
end

function d = load_cube_phiscan(cubedir, runs, scan_var_name, info, geometry)
    d = Cube(cubedir, runs, scan_var_name, info); % instantiate and construct Cube object
    d.norm_i0(); % divide by i0 values (accounting for XFEL intensity fluctuation)
    d.add_geometry(geometry); % configure the geometry struct for hkl mapping
end