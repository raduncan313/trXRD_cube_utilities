function [P, tau_p, attypes_p, atnums_p, ifc, elements] = readVASPfiles(poscar, sposcar, force_file)

ry_to_thz           = 3289.842;
amu_to_ry           = 911.4442;
bohr_angstr         = 0.52917720859;
ry_to_ev            = 13.6057;

[P, tau_p, attypes_p, atnums_p] = readPOSCAR(poscar);
[A_s, tau_s, atnums_s] = readSPOSCAR(sposcar);
FC = readFORCE(force_file);

%% convert to cartesian coordinates [Angstroms]
tau_s = A_s*tau_s;                  % atomic coord in the supercell [Angstroms]
tau_p = P*tau_p;                    % atomic coord in the unit cell [Angstroms]

nat_p = size(tau_p, 2);
nat_s = size(tau_s,2);
nuc = nat_s/nat_p;                  % number of unit cells

ifc=struct();
ifc.tau_diffs   = {};
ifc.num_uc      = nuc;
ifc.num_at      = nat_p;
ifc.ifc         = zeros(3,3,nat_p,nat_p,nuc);

%%
% find n2
centroid_s = mean(tau_s,2);
Rcen = centroid_s;
% coords= A_s\centroid_s;
% coords = coords - round(coords);
% Rcen = A_s*coords;
[dis, n2] = min( sum((Rcen - tau_s).^2, 1));

%%
% n2 = floor(n2/nuc);
% n2 =2;
n2 = floor(nuc/2);
% n2 = 2;

kk=1;
for n1=1:ifc.num_uc
    for nb=1:ifc.num_at
        for na=1:ifc.num_at
            RR = tau_s(:,n1+nuc*(na-1)) - tau_s(:,n2+nuc*(nb-1));
            [shortvecs, ~] = findinwscell(RR, A_s, 1e-6);
            for ii = 1:size(shortvecs,2)
                shortvecs(:,ii) = A_s*shortvecs(:,ii);
            end
            ifc.tau_diffs{kk} = shortvecs;   %  push!(vaspifc.tau_diffs, shortvecs)
            ifc.ifc(:,:,na,nb,n1) = FC(:,:,n1+nuc*(na-1), n2+nuc*(nb-1));
            kk = kk+1;
        end
    end
end

alltaudiffs = [];
multis = [];
for ii=1:length(ifc.tau_diffs)
    ln = size(ifc.tau_diffs{ii}, 2);
    multis = [multis, ln];
%     for k=1:ln
    alltaudiffs = [alltaudiffs, ifc.tau_diffs{ii}];
%     end
end
%
elements = [];
for ii=1:length(attypes_p)
    elements = [elements; repmat(attypes_p(ii), atnums_p(ii),1)];
end

% # to convert to Ry/bohr^2 multiply by bohr_angstr^2/ry_to_ev
ifc.ifc = ifc.ifc*bohr_angstr^2/ry_to_ev;

ifc.alltaudiffs = alltaudiffs;
ifc.multiplicities = multis;
