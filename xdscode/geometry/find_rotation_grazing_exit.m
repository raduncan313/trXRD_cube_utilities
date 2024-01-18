function Rtot = find_rotation_grazing_exit(geometry, sample_normal, hkl, beta,x0)
%FIND_ANGLES_GRAZING_EXIT finds sample and outgoing wavevector angles.
%     For a given set of miller indices `hkl`, the function finds the Huber
%     diffractometer angles, and the angles of the outgoing wavevector
%     such that the Bragg reflection is satisfied.
%
%     [phi, theta, chi, del,mu] = FIND_ANGLES_GRAZING_EXIT(geometry, hkl, beta)
%     uses the geometry structure 'geometry' to pass parameters (lattice
%     constants, x-ray energy, etc). HKL is the miller indices (in the
%     lattice defined in geometry, and `beta` is the desired
%     GRAZING EXIT ANGLE. The output, phi, theta, and chi, are the Huber  
%     diffractometer angles that satisfy the Bragg condition at such incidence angle.
%
%
%%
lambda0=geometry.lambda0;
bg = geometry.primvects';
hkl=hkl(:);
SamRot = geometry.SamRot;

%% vectors in reciprocal space
Ghkl = SamRot*bg*hkl;
k0 = 1/lambda0*[1,0,0]';

%% Bragg angle
thetaB = asind(lambda0*norm(Ghkl)/2);

%% scattering angle
tth = 2*thetaB;

%% angle between Ghkl and the incident direction
% gam1 = acosd(nsample*Ghkl./norm(Ghkl));
gam1 = acosd(k0'*Ghkl./norm(Ghkl)/norm(k0));

%% this makes the Bragg condition for a tentative intial configuration
% rotate around an axis perp to `Ghkl` and `k0`
newaxis=cross(k0, Ghkl);

% find the rotation matrix for the Bragg condition
MBragg=rotationmat3D(-gam1+90+thetaB, newaxis);

% and the rotated vectors are:
rotGhkl = MBragg*Ghkl;
rotNsample = MBragg*SamRot*sample_normal(:);

% the scattered vector is:
kp = rotGhkl+k0;

% whose norm should be the same as norm(n0)
norm(kp) - norm(k0);

%% find correct rotation to satisfy exit angle requirement
% Important, any rotation around rotGhkl maintains the Bragg condition. 
% We need to find a rotation around rotGhkl that satisfies
% kp'*rotNsample = sin(beta). 
% Note that we use a ficticious angle, not a typical lab angle. We will map
% it to Huber angles later once we have the full rotation matrix.


%% zeros of this function are the desired condition
function y=tmpfunc(ph)
%     y = norm(k0 + rot_matrix(ph, beta, chi)*SamRot*Ghkl) - norm(k0);
    MM = rotationmat3D(ph, rotGhkl);
    rotn = MM*rotNsample(:);
    y = sind(beta) - kp(:)'*rotn(:)/norm(kp);
end

%% find phi that satisfies the Bragg condition
opts=optimset('Display','off');
[phi0, f0] = fzero(@tmpfunc, x0, opts);
% phi0 
warning on
if isnan(phi0)
    warning('phi is NaN, was not able to satisfy the exit condition')
end
warning off

%% get the final rotation
if ~isnan(phi0)
    Rtot = rotationmat3D(phi0, rotGhkl)*MBragg;
else
    Rtot = MBragg;
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