%%
addpath(genpath('../../'))
% cd '/Users/mtrigo/Documents/THz-STO-KTO/STO/'

ry_to_thz=3289.842;
amu_to_ry=911.4442;

load sto_cubic.mat % load geometry and ifc structures

poscar = 'STO_300K/POSCAR';
sposcar = 'STO_300K/SPOSCAR';
force_file='STO_300K/FORCE_CONSTANTS';


% define lattice
at = geometry.realvecs';
bg = geometry.primvects';
tau = geometry.basis';

%% read VASP files, these data should be in sto_cubic.mat already, but for testing:

[P, tau_p, attypes_p, atnums_p] = readPOSCAR(poscar);
[A_s, tau_s, atnums_s] = readSPOSCAR(sposcar);
FC = readFORCE(force_file);
tau_s = A_s*tau_s;                  % atomic coord in the supercell [Angstroms]
tau_p = P*tau_p;                    % atomic coord in the unit cell [Angstroms]
[P, tau_p, attypes_p, atnums_p, ifc,elements] = readVASPfiles(poscar, sposcar, force_file);
geometry.element = elements;
geometry.realvecs = P;
geometry.primvects = inv(P')';
geometry.basis = tau_p'*P;

%% plot the coordinates of what we read
figure(4); clf

plotunitcell(P);
hold on
plotunitcell(A_s);
hold all
scatter3(tau_p(1,:),tau_p(2,:), tau_p(3,:),150,'filled','r')
scatter3(tau_s(1,:),tau_s(2,:), tau_s(3,:),'filled','b')
% scatter3(centroid_s(1,:),centroid_s(2,:), centroid_s(3,:),'filled','g')
% scatter3(Rcen(1,:),Rcen(2,:), Rcen(3,:),200,'k')

axis vis3d
rotate3d on
view(20, 45)
% 
%% set up a q-path along a high symmetry path
at = geometry.realvecs';
bg = geometry.primvects';

hklpts = [ 0 0 0; 0 1 0; 1 1 0; 0 0 0; 1 1 1; 0 1 0]'/2;
kxyzpts = bg*hklpts;
kxyzdiff = zeros(3,size(kxyzpts, 2)-1);

for ii = 1:size(kxyzdiff, 2)
    kxyzdiff(:,ii) = kxyzpts(:,ii+1) - kxyzpts(:, ii);
end

npts = 100;

qpath = zeros(3,npts*size(kxyzdiff, 2)+1);
flngth = linspace(0,1-1/npts,npts);
%
for ii = 1:size(kxyzdiff, 2)
    qpath(:,npts*(ii-1)+1:npts*ii) = ...
        bsxfun(@plus, kron(flngth, kxyzdiff(:,ii)), kxyzpts(:,ii));
end
qpath(:,end) = kxyzpts(:,end);
tpts = 1:npts:size(qpath,2);
ptnames = {'\Gamma', 'X', 'M', '\Gamma',  'R', 'X'};

%% compute and plot the phonon dispersion
DD2 = dynmat2(qpath', ifc, geometry.masses*amu_to_ry);
[ww, vv] = phonon_freqs_vecs(DD2);

figure(2); clf
plot(ry_to_thz*ww','r')
set(gca, 'XTick', tpts)
set(gca, 'XTickLabel', ptnames)
axis tight

%% alternative way to compute the TDS intensities, frequencies and vectors
[intens, ww, vv] = tds_intensity_onq(qpath', geometry, ifc);

figure(3); clf
plot(ww','r')
set(gca, 'XTick', tpts)
set(gca, 'XTickLabel', ptnames)
axis tight


%% plot STO diffuse scattering in a reasonable experimental geometry

geometry.detector.det_dist=100;
geometry.detector.det_pixels_horz=2*1024;
geometry.detector.det_pixels_vert=2*1024;
geometry.detector.det_size_horz = 0.070*geometry.detector.det_pixels_horz;
geometry.detector.det_size_vert = 0.070*geometry.detector.det_pixels_vert;
geometry.imageNy=128;
geometry.imageNz=128;

detRot = rotationmat3D(-90, [0,1,0]);
geometry.detector.detRot = detRot;
geometry.rot_matrix = @huber_matrix;

geometry.SamRot = huber_matrix(0, 0, 0);
geometry.phi = 40;
geometry.theta= 1;

masses = geometry.masses*amu_to_ry;
nat=ifc.num_at;
nuc=ifc.num_uc;
bg = geometry.primvects';
at = geometry.realvecs';

% ************************************************************************
% calculate only the TDS intensity on the Ewald sphere, returns the frequencies 
% and eigenvectors for free
%
[intens, ww, vv] = tds_intensity_grid(geometry, ifc);

% plot the intensity
cc = reshape(sum(intens,1), geometry.imageNz, geometry.imageNy);
figure(4); clf
plotdata(geometry, (cc));
colorbar

%% ************************************************************************
% calculate dynmats, eigenfrequencies and vectors on the Ewald sphere
%

Q=det_kspace_sphere(geometry);
Qflat=reshape(Q,[],3)*geometry.rot_matrix(geometry.phi,geometry.theta, geometry.chi)*geometry.SamRot;
nq=size(Qflat,1);

% tic
DD2 = dynmat2(Qflat, ifc, geometry.masses*amu_to_ry);
[ww, vv] = phonon_freqs_vecs(DD2);

% compute TDS intensity on 2D grid
I11=tdsI1(geometry, Qflat, ww, vv);
cc = reshape(sum(I11,1), geometry.imageNz, geometry.imageNy);
figure(3); clf
plotdata(geometry, (cc));
