%% setup
addpath(genpath('../../'))

ry_to_thz=3289.842;
amu_to_ry=911.4442;

% has geometry structure and IFCs 
s = load('snse.mat');
geometry = s.geometry;

% at=geometry.realvecs';
% geometry.basis = at*geometry.basis;

poscar = 'POSCAR';
sposcar = 'SPOSCAR';
force_file='FORCE_CONSTANTS';
% 

%% read VASP files

[P, tau_p, attypes_p, atnums_p] = readPOSCAR(poscar);
[A_s, tau_s, atnums_s] = readSPOSCAR(sposcar);
FC = readFORCE(force_file);
tau_s = A_s*tau_s;                  % atomic coord in the supercell [Angstroms]
tau_p = P*tau_p;                    % atomic coord in the unit cell [Angstroms]
[P, tau_p, attypes_p, atnums_p, ifc] = readVASPfiles(poscar, sposcar, force_file);

%% set up a q-path along high symmetry path
% TODO: fix the labels and q-points
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

figure(2); clf
plot(ry_to_thz*ww')
set(gca, 'XTick', tpts)
set(gca, 'XTickLabel', ptnames)
ylabel('Frequency (THz)')
xlabel('wavevector')

