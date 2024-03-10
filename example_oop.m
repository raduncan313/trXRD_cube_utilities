clear all
close all

addpath(genpath('xdscode'));
cubedir = 'C:/Users/radun/OneDrive - Stanford/Desktop/TSI_SACLA_June2023/cubes_post/';
t0 = 4; % Nominal time-zero for this dataset.
runs = {1292991};
info = 'testing';
% d = Cube(cubedir, runs, info);
d = Cube('./', {129}, 'testing');
d.norm_i0();
f = d.plotdatacube('on', [0 100], []);
% f = d.plot_lineout();

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