function Rtot = find_rotation_grazing_incidence(geometry, hkl, sample_normal, alpha, x0)
%find_rotation_grazing_incidence finds sample and outgoing wavevector angles.
% Finds the rotation matrix to satify the Bragg condition at an angle of
% incidence with respect to the sample surface (specified in
% `sample_normal`) given by `alpha`.
%
%%
lambda0=geometry.lambda0;
bg = geometry.primvects';
hkl=hkl(:);
SamRot = geometry.SamRot;

%% vectors in reciprocal space
Ghkl = SamRot*bg*hkl;
k0 = 1/lambda0*[1,0,0]';

%% vectors in real space
sample_normal = SamRot*sample_normal(:);

%% angle between incident direction and sample normal
alph1= acosd(sample_normal(:)'*k0(:)./norm(k0)/norm(sample_normal));

%% rotate around this axis to satisfy the incidence condition
newaxis=cross(k0, sample_normal);
Minc=rotationmat3D(-alph1+90+alpha, newaxis);

% and the rotated vectors in the grazing configuration are:
rotGhkl = Minc*Ghkl;
rotNsample = Minc*sample_normal(:);

%% Find rotation for Bragg condition

%% zeros of this function are the desired condition
function y=tmpfunc(ph)
%     y = norm(k0 + rot_matrix(ph, alpha, chi)*SamRot*Ghkl) - norm(k0);
    MBragg = rotationmat3D(ph, rotNsample);
    y = norm(k0 + MBragg*rotGhkl) - norm(k0);
end

%% find phi that satisfies the Bragg condition
opts=optimset('Display','off');
[phi0, f0] = fzero(@tmpfunc, x0, opts);
% phi0 
warning on
if isnan(phi0)
    warning('phi is NaN, was not able to satisfy the Bragg condition')
end
warning off

%% get the final rotation
if ~isnan(phi0)
    Rtot = rotationmat3D(phi0, rotNsample)*Minc;
else
    Rtot = Minc;
end

%% if the normal and k0 are on the same side (Laue) try to find a new solution 
% at 180 deg.
% if k0(:)'*(Rtot*rotNsample(:)) > 0
%     [phi0, f0] = fzero(@tmpfunc, 90, opts);
%     Rtot = rotationmat3D(phi0, rotGhkl)*MBragg;
% end

%
% and the final rotated vectors are:
% finalRotGhkl    = Rtot*Ghkl;
% finalRotNsample = nsample*Rtot';
% finalkp         = k0 + finalRotGhkl;


end