clear all
close all

% This is an example script for how to use the various functions in this
% repository. Ask me (Ryan) for the data needed to run this.

cubedir = 'C:/Users/radun/OneDrive - Stanford/Desktop/TSI_SACLA_June2023/cubes_post/';
t0 = 4; % Nominal time-zero for this dataset.

r1l_4m31 = {1292991};
d1 = load_and_preprocess(cubedir, r1l_4m31, t0, '4m31 1 mJpcm2 long ROIs'); % struct for saving ROIs
d2 = load_and_preprocess(cubedir, r1l_4m31, t0, '4m31 1 mJpcm2 long threshs'); % struct for saving thresholded signal

%% Plot ROIs and save signals
[f1, d1] = plot_rois(d1, 2)
write_sigs_to_csv(d1, './sigs')

%% Plot sum of pixels with total intensity over a given threshold (relative to the max pixel intensity) and save signals
[f2, d2] = thresh_and_plot(d2, 0.05)
write_sigs_to_csv(d2, './sigs')

%% Plot a lineout over the detector
[f3, d1] = plot_lineout(d1);

%% Helper functions for preprocessing

function d = load_and_preprocess(cubedir, runs, t0, info)
    d = read_cubes(cubedir, runs);
    d = norm_i0(d);
    d = orient_SACLA(d);
    d.on.scan_var = d.on.scan_var - t0;
    d.off.scan_var = d.off.scan_var - t0;
    d.info = info;
end

function d = orient_SACLA(d)
    d.on.imgs = permute(d.on.imgs, [2,1,3]);
    d.off.imgs = permute(d.off.imgs, [2,1,3]);
end