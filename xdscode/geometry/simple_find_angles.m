function [phi, theta, chi, mu, del] = simple_find_angles(geometry, hkl, theta, theta_start)
% function [phi, theta, chi, del, mu] = find_angles_outgoingk(geometry, hkl, theta)
%FIND_ANGLES_OUTGOINGK finds Outgoing photon wavevector angles.
%     For a given set of miller indices hkl, the function finds the Huber
%     diffractometer angles such that the Bragg reflection is satisfied.
%
%     [phi, theta, chi] = FIND_ANGLES_OUTGOINGK(geometry, hkl, theta)
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

% should check this is not empty:
rot_matrix = geometry.rot_matrix;

%% vectors in reciprocal space
Ghkl = bg*hkl;
k0 = 1/lambda0*[1,0,0];

%% Bragg angle
thetaB = asind(lambda0*norm(Ghkl)/2);

%% scattering angle
tth = 2*thetaB;

%% angle between Ghkl and the sample normal
gam1 = acosd([0, 1, 0]*Ghkl./norm(Ghkl));

%% 
% Rsample = rotationmat3D(90, [1,0,0]);
% Ghkl = Rsample*Ghkl(:);   % point b-axis along z-axis
% Ghkl    % in this geometry b-axis points along y (horizontal)
Ghkl = Ghkl(:);

%% the angle between the surface and the incident beam
alph1 = thetaB - gam1;
chi=0;

%% zeros of this function are the Bragg condition
function y=tmpfunc(ph)
    k0=k0(:);
    y = norm(k0 + rot_matrix(ph, theta, chi)*SamRot*Ghkl) - norm(k0);
end

%% find phi that satisfies the Bragg condition
opts=optimset('Display','off');
[phi0, f0] = fzero(@tmpfunc, theta_start, opts);
% phi0 

%% SPECIAL CASE when the above cannot be satisfied (e.g. symmetric reflection)
% if we cannot find a solution, display the Bragg angle and incidence angle
% for "in-plane" scattering
if isnan(phi0)
    disp(['Cannot satisfy the Bragg condition, ',...
        'assuming symmetric Bragg...']);
    
    % set Theta=alph1
    theta = alph1;
%     Ghkl_parallel = Ghkl(:).*[1,1,0]';
    phi = 180-180/pi*atan2(Ghkl(2),Ghkl(1));
    chi = 0;
    % get Rtot from these
    Rtot = rot_matrix(phi, theta, chi);
    % 
else
%% rotation matrix for such phi
    Rtot = rot_matrix(phi0, theta, chi);
end

%% check that Bragg condition is satisfied with this

%rotated Ghkl
rotGhkl = Rtot*SamRot*Ghkl(:);

% k-prime defined as
kp=rotGhkl+k0(:);

% this should be zero
norm(k0)-norm(kp);

%% calculate the outgoing wavevector angles

% angle in the horizontal plane
% kp1=kp-([0,0,1]*kp)*kp/norm(kp)

% angle in the vertical plane
% kp2=kp-([0,1,0]*kp)*kp/norm(kp)

%polar coordinates
mu = atan2d(kp(2), kp(1));

% mu = atand(kp1(2)/kp1(1));
del = 90.0 - acosd(kp(3)*lambda0);
% mu = acosd(kp1'*k0(:)/norm(kp1)/norm(k0));

%% final full rotation matrix
% Rtot = rotationmat3D(roll_angle, k0)*huber_rotation(phi0, theta, chi);
% Rtot = rotGhkl;

%% Get the Huber angles from the rotation matrix above
% 
phi     = 180/pi*atan2(Rtot(2,1),Rtot(2,2));
theta   = 180/pi*atan2(-Rtot(1,3),Rtot(3,3));
chi     = asind(-Rtot(2,3));

%% print diagnostics and output

% fprintf(['\n\n\nScattering vector in lab coordinates:',... 
% '\nG(h,k,l) = [%f, %f, %f]\n'], Rtot*Ghkl(:))

% fhkl2 = abs(Fhkl(geometry, hkl, 0*geometry.basis))^2;

fprintf(['HKL\t\tThetaB\t|Phi\tTheta\tMu\tDelta|\n'])
fprintf(['%2.2f %2.2f %2.2f\t%2.3f\t'], hkl, thetaB)
fprintf(['%2.3f  %2.3f  %2.3f  %2.3f\n\n'], phi, theta, mu, del)
% n1=Rsample*[0,1,0]';
% alph1=acosd(k0'*Rtot*n1/norm(k0))-90;

end
