function [phi, theta, chi, tth, fhkl2] = find_angles_planar_geometry(geometry, hkl, theta, init_phi)
%FIND_ANGLES_PLANAR_GEOMETRY finds the Huber angles to meet the diffraction condition.
%     For a given set of miller indices hkl, the function finds the Huber
%     diffractometer angles such that the Bragg reflection is satisfied and
%     the reflected beam is in the horizontal plane.
%
%     [phi, theta, chi, tth, fhkl2] = FIND_ANGLES_PLANAR_GEOMETRY(geometry, hkl, theta)
%     uses the geometry structure 'geometry' to pass parameters (lattice
%     constants, x-ray energy, etc). HKL is the miller indices (in the
%     lattice defined in geometry, and the input theta is the desired
%     INCIDENCE ANGLE (not the Bragg angle). The output, phi, theta, and
%     chi, are the Huber diffractometer angles that satisfy the Bragg
%     condition at such incidence angle.
%
%%

lambda0=geometry.lambda0;
bg = geometry.primvects';
hkl=hkl(:);
SamRot = geometry.SamRot;

% we assume that initially there is no Chi rotation:
chi=0;

%% vectors in reciprocal space
Ghkl = bg*hkl;
k0 = 1/lambda0*[1,0,0];

%% Bragg angle
thetaB = asind(lambda0*norm(Ghkl)/2);

%% scattering angle
tth = 2*thetaB;

%% angle between Ghkl and the sample normal
gam1 = acosd([0, 1, 0]*Ghkl./norm(Ghkl));

%% the angle between the surface and the incident beam
alph1 = thetaB - gam1;

%% zeros of this function are the Bragg condition
function y=tmpfunc(ph)
    k0=k0(:);
    y = norm(k0 + huber_matrix(ph, theta, chi)*SamRot*Ghkl) - norm(k0);
end

%% find phi that satisfies the Bragg condition
opts=optimset('Display','off');
[phi0, f0] = fzero(@tmpfunc, init_phi, opts);
% phi0 

%% SPECIAL CASE when the above cannot be satisfied (e.g. symmetric reflection)
% if we cannot find a solution, display the Bragg angle and incidence angle
% for "in-plane" scattering
if isnan(phi0)
    warning on
    warning(['Cannot satisfy the Bragg condition, ',...
        'assuming symmetric Bragg...'])
    
    % set Theta=alph1
    theta = alph1;
%     Ghkl_parallel = Ghkl(:).*[1,1,0]';
    phi = 180-180/pi*atan2(Ghkl(2),Ghkl(1));
    chi = 0;
    % get Rtot from these
    RHuber = huber_matrix(phi, theta, chi);
    % 
else
%% rotation matrix for such phi
    RHuber = huber_matrix(phi0, theta, chi);
end

%% check that Bragg condition is satisfied with this

%rotated Ghkl
rotGhkl = RHuber*SamRot*Ghkl(:);

% k-prime defined as
kp=rotGhkl+k0(:);

% this should be zero
norm(k0)-norm(kp);

%% now we want to bring the reflection to the horizontal plane
% 
% this is the roll angle around the incident k-vector
roll_angle = -180/pi*atan2(rotGhkl(3), rotGhkl(2));

%% final full rotation matrix
Rtot = rotationmat3D(roll_angle, [1,0,0])*huber_matrix(phi0, theta, chi);

%% Get the Huber angles from the rotation matrix above
phi     = 180/pi*atan2(Rtot(2,1),Rtot(2,2));
theta   = 180/pi*atan2(-Rtot(1,3),Rtot(3,3));
chi     = asind(-Rtot(2,3));

%% check that the scattering vector is in the plane
fhkl2 = abs(Fhkl(geometry, hkl))^2;
% fprintf('**************************************\n')
% fprintf('Output:\n')
fprintf('HKL        ThetaB    |Phi    Theta    Chi|    |Fhkl|^2\n')
fprintf('%i %i %i\t    %2.3f   %f   ', hkl, thetaB)
fprintf('%2.3f   %2.3f   %2.3f   %2.3f\n', phi, theta, chi, fhkl2)
% fprintf('Scattering vector in lab coordinates: (should be ~ [X, Y, 0])\n')
% fprintf('G(h,k,l) = [%f, %f, %f]\n', Rtot*SamRot*Ghkl(:))
% fprintf('**************************************\n')

n1=SamRot\[0,0,1]';
alph1=acosd(k0'*Rtot*n1/norm(k0))-90;

end