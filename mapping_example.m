clear all
close all

addpath(genpath('xdscode'))
load r129.mat
[P, tau_p, attypes_p, atnums_p] = readPOSCAR('TaSe42I_sm_isp_sd_1502042_springermaterials.vasp');

tau_p = P*tau_p; % convert to cartesian
elements = [];
for ii=1:length(attypes_p)
    elements = [elements; repmat(attypes_p(ii), atnums_p(ii),1)];
end
geometry.element = elements;

at = P;
bg = inv(at);
geometry.realvecs = at';
geometry.primvects = bg';
geometry.basis=tau_p';
geometry.lambda0 = 12.398/9.55;

k0 = 1/geometry.lambda0*[1,0,0];

for ii = 1:length(geometry.element)
    el = geometry.element{ii};
    geometry.form_factors{ii} = make_xrayscatt_factor(el, 12398./geometry.lambda0);
end

geometry.detector.det_dist = 600;
geometry.detector.det_pixels_horz=1030;
geometry.detector.det_pixels_vert=1064;
geometry.detector.det_size_horz = 0.075*geometry.detector.det_pixels_horz;
geometry.detector.det_size_vert = 0.075*geometry.detector.det_pixels_vert;
geometry.imageNy=round(geometry.detector.det_pixels_horz/10);
geometry.imageNz=round(geometry.detector.det_pixels_vert/10);
geometry.beam_center = [580, 534];

geometry.theta = -0.5;
geometry.phi = -22.2;
geometry.chi = 90;
geometry.mu = 0;
geometry.nu = -24.57;
geometry.delta = 20.6;


sample_normal=[1,1,0];
sample_inplane=[0,0,1];
geometry.SamRot = rotationmat3D(45,[0,0,1]);
geometry.rot_matrix = @(a,b,c,d) huber_sixcircle_matrix_lw89(a,b,c,d);

ROI_lim_h = [300 550];
ROI_lim_v = [100 500];
ROI_lims = {ROI_lim_h, ROI_lim_v};
scan_range = linspace(-4, 4, length(cube(1,1,:)));
scan_angle = 'phi';
threshs = [40 1000];

hkl_scatter = map_huber_anglescan_to_hkl(cube, ROI_lims, scan_angle, scan_range, threshs, geometry);
scatter3(hkl_scatter(:,1),hkl_scatter(:,2),hkl_scatter(:,3),10*ones(length(hkl_scatter(:,1)),1),hkl_scatter(:,5),'filled')
xlabel('h')
ylabel('k')
zlabel('l')
axis vis3d


