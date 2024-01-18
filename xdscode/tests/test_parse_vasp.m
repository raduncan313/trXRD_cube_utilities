%%
addpath(genpath('~/Dropbox/xdscode/'))
cd ~/Dropbox/xdscode/tests/

ry_to_thz           = 3289.842;
amu_to_ry           = 911.4442;
bohr_angstr         = 0.52917720859;
ry_to_ev            = 13.6057;

s = load('vo2m1.mat','geometry');
% s = load('vo2rutile.mat','geometry');
geometry = s.geometry;

% poscar = 'mp-4651-20180417/POSCAR';
% sposcar = 'mp-4651-20180417/SPOSCAR';
% force_file='mp-4651-20180417/FORCE_CONSTANTS';
% 
poscar = '../examples/VO2/M1/POSCAR';
sposcar = '../examples/VO2/M1/SPOSCAR';
force_file='../examples/VO2/M1/FORCE_CONSTANTS';
% 
% poscar = '../examples/VO2/rutile/POSCAR';
% sposcar = '../examples/VO2/rutile/SPOSCAR';
% force_file='../examples/VO2/rutile/FORCE_CONSTANTS';

[P, tau_p, attypes_p, atnums_p] = readPOSCAR(poscar);
[A_s, tau_s, atnums_s] = readSPOSCAR(sposcar);
FC = readFORCE(force_file);

%% convert to cartesian coordinates [Angstroms]
tau_s = A_s*tau_s;                  % atomic coord in the supercell [Angstroms]
tau_p = P*tau_p;                    % atomic coord in the unit cell [Angstroms]

nat_p = size(tau_p, 2);
nat_s = size(tau_s,2);
nuc = size(tau_s,2)/nat_p;          % number of unit cells

clear ifc
% tau_diffs = zeros(3,3,nat_p,nat_p,nuc);
ifc=struct();
ifc.tau_diffs   = {};
ifc.num_uc      = nuc;
ifc.num_at      = nat_p;
ifc.ifc         = zeros(3,3,nat_p,nat_p,nuc);

%% this does not produce what we want
% find n2
% centroid_s = mean(tau_s,2);
% Rcen = centroid_s;
% coords= A_s\centroid_s;
% coords = coords - round(coords);
% Rcen = A_s*coords;
% [dis, n2] = min( sum((Rcen - tau_s).^2, 1));

%%
n2 = floor(nuc/2); % a good guess
kk=1;
for n1=1:ifc.num_uc
    for nb=1:ifc.num_at
        for na=1:ifc.num_at
            RR = tau_s(:,n1+nuc*(na-1)) - tau_s(:,n2+nuc*(nb-1));
            [shortvecs, ~] = findinwscell(RR, A_s, 1e-6);
%             for ii in 1:length(shortvecs)
            for ii = 1:size(shortvecs,2)
                shortvecs(:,ii) = A_s*shortvecs(:,ii);
            end
            ifc.tau_diffs{kk} = shortvecs;   %  push!(vaspifc.tau_diffs, shortvecs)
            ifc.ifc(:,:,na,nb,n1) = FC(:,:,n1+nuc*(na-1), n2+nuc*(nb-1));
            kk = kk+1;
        end
    end
end
% # to convert to Ry/bohr^2 multiply by bohr_angstr^2/ry_to_ev
ifc.ifc = ifc.ifc*bohr_angstr^2/ry_to_ev;

alltaudiffs = [];
multis = [];
for ii=1:length(ifc.tau_diffs)
    ln = size(ifc.tau_diffs{ii}, 2);
    multis = [multis, ln];
%     for k=1:ln
    alltaudiffs = [alltaudiffs, ifc.tau_diffs{ii}];
%     end
end

ifc.alltaudiffs = alltaudiffs;
ifc.multiplicities = multis;

%%
elements = [];
for ii=1:length(attypes_p)
    elements = [elements; repmat(attypes_p{ii}, atnums_p(ii),1)];
end


%% cells above became a function:
[P, tau_p, attypes_p, atnums_p, ifc, elements] = readVASPfiles(poscar, sposcar, force_file);

%% set up a q-path along high symmetry path
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

%% compute diffuse scattering
geometry.rot_matrix=@huber_matrix;
% the next few cells define the geometry parameters and orientation for the
% various possible VO2 twins
% theta=.5;
% phi=5.7-7+22.;
geometry.phi=1;

nqz = 64;
nqy = 64;
geometry.imageNz = nqz;
geometry.imageNy = nqy;

%
[intens, ~, ~] = tds_intensity_grid(geometry, ifc);
cc = reshape(sum(intens,1), geometry.imageNz, geometry.imageNy);
figure(1);clf
plotdata(geometry, cc, [0,1]);
caxis auto

%% plot the coordinates of what we read
figure(4); clf

plotunitcell(P);
hold on
plotunitcell(A_s);
hold all
scatter3(tau_p(1,:),tau_p(2,:), tau_p(3,:),150,'filled','r')
scatter3(tau_s(1,:),tau_s(2,:), tau_s(3,:),'filled','b')
scatter3(centroid_s(1,:),centroid_s(2,:), centroid_s(3,:),'filled','g')
scatter3(Rcen(1,:),Rcen(2,:), Rcen(3,:),200,'k')

axis vis3d
rotate3d on
view(20, 45)
% 

%%
% masses = geometry.masses*amu_to_ry;
% nat=ifc.num_at;
% nuc=ifc.num_uc;
% bg=geometry.primvects';
% 
% % flatten the tau_diffs cache of vectors         ----generic----
% tic
% alltaudiffs = [];
% multis = [];
% for ii=1:length(ifc.tau_diffs)
%     ln = length(ifc.tau_diffs{ii});
%     multis = [multis, ln];
%     for k=1:ln
%         alltaudiffs = [alltaudiffs, ifc.tau_diffs{ii}{k}];
%     end
% end
% toc
% ifc.alltaudiffs = alltaudiffs;
% ifc.multiplicities = multis;

