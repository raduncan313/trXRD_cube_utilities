clear all
close all

cubedir = 'C:/Users/radun/OneDrive - Stanford/Desktop/TSI_SACLA_June2023/cubes_post/';
t0 = 4; % Nominal time-zero for this dataset.
runs = {1292991};
info = 'testing';
d1 = Cube(cubedir, runs, info);
d1.norm_i0();
d1.thresh_and_plot(0.05);
d1.write_sigs_to_csv('sigs/')

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