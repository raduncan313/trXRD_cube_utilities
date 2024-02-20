clear all
close all

% This is an example script for how to use the various functions in this
% repository. Ask me (Ryan) for the data needed to run this.

cubedir = 'C:/Users/radun/OneDrive - Stanford/Desktop/TSI_SACLA_June2023/cubes_post/';
csvdir = 'sigs/';
th = 0.05; % Pixel intensity treshold (relative to max pixel intensity) for defining the ROI
t0 = 4.3; % Nominal time-zero (needs a bit of tweaking on a run-by-run basis)
cmap = lines; % colormap used later for plotting

r1l_4m31 = {1292991};
d1 = load_and_preprocess(cubedir, r1l_4m31, th, t0 - 0.3, '4m31 1 mJpcm2 long ROIs', 1); % struct for saving ROIs
d2 = load_and_preprocess(cubedir, r1l_4m31, th, t0 - 0.3, '4m31 1 mJpcm2 long threshs', 1); % struct for saving thresholded signal

%% Plot ROIs and save signals
[f1, d1] = plot_rois(d1, 2)
write_sigs_to_csv(d1, './sigs')

%% Plot sum of pixels with total intensity over a given threshold (relative to the max pixel intensity) and save signals
[f2, d2] = thresh_and_plot(d2, 0.05)
write_sigs_to_csv(d2, './sigs')

%% Plot a lineout over the detector
[f3, d1] = plot_lineout(d1);

%% Helper functions for preprocessing

function d = load_and_preprocess(cubedir, runs, th, t0, info, col_ind)
    d = read_cubes(cubedir, runs);
%     d = rebin_cube(d, 1);
    d = norm_i0(d);
    d = orient_SACLA(d);
    % d = thresh_mask(d, th);
    d.on.scan_var = d.on.scan_var - t0;
    d.off.scan_var = d.off.scan_var - t0;
    d.info = info;
    % d.col_ind = col_ind;
end

function d = orient_SACLA(d)
    d.on.imgs = permute(d.on.imgs, [2,1,3]);
    d.off.imgs = permute(d.off.imgs, [2,1,3]);
end