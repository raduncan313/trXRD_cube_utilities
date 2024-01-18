%% setup
addpath(genpath('../../'));

ry_to_thz=3289.842;
amu_to_ry=911.4442;

load ge.mat  % has geometry structure and IFCs 
% load('si.mat')

%% tweak the geometry parameters
%

% Detector size and coordinates:
geometry.detector.det_dist=195;
geometry.detector.det_pixels_horz=1750;
geometry.detector.det_pixels_vert=1750;
geometry.detector.det_size_horz = 198;
geometry.detector.det_size_vert = 198;
geometry.detector.detRot        = eye(3);
geometry.beam_center = [1398.9 2090];


% set the sample orientation using the Huber diffractometer angles:
geometry.SamRot = huber_matrix(0, 0, 0);
geometry.phi = 22.5;
geometry.theta= 0.2;

masses = geometry.masses*amu_to_ry;
nat=ifc.num_at;
nuc=ifc.num_uc;
bg=geometry.primvects';

% setup for internal representation of force constants
alltaudiffs = [];
multis = [];
for ii=1:length(ifc.tau_diffs)
    ln = length(ifc.tau_diffs{ii});
    multis = [multis, ln];
    for k=1:ln
        alltaudiffs = [alltaudiffs, ifc.tau_diffs{ii}{k}];
    end
end
ifc.alltaudiffs = alltaudiffs;
ifc.multiplicities = multis;

%% ************************************************************************
% calculate TDS pattern for given sample orientation 
%
theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;
Q=det_kspace_proj(geometry);
Qflat=reshape(Q,[],3)*Rot;
nq=size(Qflat,1);

% tic
DD2 = dynmat2(Qflat, ifc, masses);
[ww, vv] = phonon_freqs_vecs(DD2);

% compute TDS intensity on 2D grid
I11=tdsI1(geometry, Qflat, ww, vv);
cc = reshape(sum(I11,1), geometry.imageNz, geometry.imageNy);
figure(3); clf
plotdata(geometry, log10(cc));
% colorbar

%% plots the frequencies on 2D grid
ww=reshape(ww, 3*nat, sqrt(nq),sqrt(nq));
freqs=ry_to_thz*ww;
figure(5); clf
plotdata(geometry, squeeze(freqs(1,:,:)), [0,3.]);
colorbar

%% plot lineouts with `qslineout`
% get the `Q`s for a lineout on the detector
[intens,qs,cx,cy]=qslineout(geometry, log10(cc),[-.2 2.5]); % uses plotdata(...)
[intens, freqs, ~] = tds_intensity_onq(qs, geometry, ifc);
figure(3);clf
subplot(2,1,1)
plot(freqs')
ylabel('Frequency (THz)')
subplot(2,1,2);
plot(qs*inv(bg)')
ylabel('Q (rlu)')
xlabel('path index')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example Data from Ge run 2015
%
load ~/Documents/xds_example/ge_data.mat

%% plots the time domain data
figure(1); clf
plotdatacube(delays, data - data(:,:,2),[], [0., 20], [0,0]);

%% take a lineout of Q-space on the experimental data
figure(4); clf
geometry.imageNy=512;
geometry.imageNz=512;
[intens,qs,cx,cy]=qslineout(geometry, (abs(data(:,:,2))),[-1 3e2]); 
[intens, freqs, ~] = tds_intensity_onq(qs, geometry, ifc);
geometry.imageNy=128;
geometry.imageNz=128;

%% produce time profiles of the qslineout chosen above
numbins = length(delays);
psps = zeros(length(cx),numbins);

% time-zero seems to be at index=25
t0id = 25; 
for ii = 1:numbins-1
    tmp = data(:,:,ii+1)-data(:,:,2);
    tmp = interp2(1:size(tmp,2),1:size(tmp,1),tmp,cx,cy);
    psps(:,ii) = tmp;
end
% figure(4);clf
% imagesc(psps); axis xy;

% FT of the indidual time-traces above:
[ws, exp_disp] = tdsfft(delays(t0id:end),psps(t0id:end,:)');

%% plot the experimental dispersion
figure(4);clf
imagesc(1:length(cx),ws,log10(abs(exp_disp))); axis xy;
caxis([-2, 0.4])
hold all
plot(1:length(cx), freqs')
