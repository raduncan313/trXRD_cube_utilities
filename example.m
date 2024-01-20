clear all
close all

cubedir = 'C:/Users/radun/OneDrive - Stanford/Desktop/TSI_SACLA_June2023/cubes_post/';
csvdir = 'sigs/';
th = 0.05; % Pixel intensity treshold (relative to max pixel intensity) for defining the ROI
t0 = 4.3; % Nominal time-zero (needs a bit of tweaking on a run-by-run basis)
cmap = lines; % colormap used later for plotting

r1s_01m9 = {1293041, 1293043, 1293044};
r3s_01m9 = {1293036, 1293039};
r7s_01m9 = {1293034};
r10s_01m9 = {1293052};
r13s_01m9 = {1293046, 1293049};

d1s_01m9 = load_and_preprocess(cubedir, r1s_01m9, th, t0 + 0.1, '01m9 1 mJpcm2 short', 1);
d3s_01m9 = load_and_preprocess(cubedir, r3s_01m9, th, t0 + 0.1, '01m9 3 mJpcm2 short', 2);
d7s_01m9 = load_and_preprocess(cubedir, r7s_01m9, th, t0 + 0.05, '01m9 7 mJpcm2 short', 3);
d10s_01m9 = load_and_preprocess(cubedir, r10s_01m9, th, t0, '01m9 10 mJpcm2 short', 5);
d13s_01m9 = load_and_preprocess(cubedir, r13s_01m9, th, t0 + 0.05, '01m9 13 mJpcm2 short', 4);

% ds = {d1s_01m9, d3s_01m9, d7s_01m9, d10s_01m9, d13s_01m9};

%%
d = d1s_01m9;

% [f, d1] = plot_rois(d, 3)
% [f, d1] = thresh_and_plot(d, 0.05)
[f, d] = plot_lineout(d);

%% Functions

function d = load_and_preprocess(cubedir, runs, th, t0, info, col_ind)
    d = read_cubes(cubedir, runs);
%     d = rebin_cube(d, 1);
    d = norm_i0(d);
    d = orient_SACLA(d);
    d = thresh_mask(d, th);
    d.on.scan_var = d.on.scan_var - t0;
    d.off.scan_var = d.off.scan_var - t0;
    d.info = info;
    d.col_ind = col_ind;
end

function d = orient_SACLA(d)
    d.on.imgs = permute(d.on.imgs, [2,1,3]);
    d.off.imgs = permute(d.off.imgs, [2,1,3]);
end