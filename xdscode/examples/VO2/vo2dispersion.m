%%
addpath(genpath('../../'))

ry_to_thz           = 3289.842;
amu_to_ry           = 911.4442;
bohr_angstr         = 0.52917720859;
ry_to_ev            = 13.6057;

s = load('vo2m1.mat','geometry');
% s = load('vo2rutile.mat','geometry');
geometry = s.geometry;

poscar = 'M1/POSCAR';
sposcar = 'M1/SPOSCAR';
force_file='M1/FORCE_CONSTANTS';
% 
%% read VASP files
[P, tau_p, attypes_p, atnums_p, ifc] = readVASPfiles(poscar, sposcar, force_file);

%% set up q-path along high symmetry path
at = geometry.realvecs';
bg = geometry.primvects';

hklpts = [0 0 1; 0 0 0; 0 1 1; 1 1 1; 1 1 0; 0 0 0; 1 0 0]'/2;
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
ptnames = {'Z', '\Gamma', 'R', 'A', 'M', '\Gamma', 'X'};

%% compute and plot the phonon dispersion
DD2 = dynmat2(qpath', ifc, geometry.masses*amu_to_ry);
[ww, vv] = phonon_freqs_vecs(DD2);

figure(1); clf
plot(ry_to_thz*ww')
set(gca, 'XTick', tpts)
set(gca, 'XTickLabel', ptnames)
ylabel('Frequency (THz)')
xlabel('wavevector')
