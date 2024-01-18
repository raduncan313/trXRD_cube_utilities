function [phi, theta, chi, mu, del] = find_angles_grazing_exit(geometry, sample_normal, hkl, beta, x0)
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

%% get the final rotation
% Rtot = rotationmat3D(phi0, rotGhkl)*MBragg;
Rtot = find_rotation_grazing_exit(geometry, sample_normal, hkl, beta, x0);
% Rtot = Rtot*SamRot;

%% vectors in reciprocal space
Ghkl = bg*hkl;
k0 = 1/lambda0*[1,0,0]';

%% Bragg angle
thetaB = asind(lambda0*norm(Ghkl)/2);

%% the final rotated vectors are:
finalRotGhkl    = Rtot*Ghkl;
finalRotNsample = Rtot*sample_normal(:);
finalkp         = k0 + finalRotGhkl;

%% check that Bragg condition is satisfied with this

%rotated Ghkl
% rotGhkl = Rtot*SamRot*Ghkl(:);

% k-prime defined as
% kp=rotGhkl+k0;

% this should be zero (Bragg)
norm(k0)-norm(finalkp);
beta - asind(finalkp(:)'/norm(finalkp)*finalRotNsample(:))  % This should be zero too

%% calculate the outgoing wavevector angles

% angle in the horizontal plane
% kp1=kp-([0,0,1]*kp)*kp/norm(kp)

% angle in the vertical plane
% kp2=kp-([0,1,0]*kp)*kp/norm(kp)

%polar coordinates
mu = atan2d(finalkp(2), finalkp(1));

% mu = atand(kp1(2)/kp1(1));
del = 90.0 - acosd(finalkp(3)*lambda0);
% del = acosd(kp1'*k0(:)/norm(kp1)/norm(k0));

%% Get the Huber angles from the rotation matrix above
% 
phi     = 180/pi*atan2(Rtot(2,1),Rtot(2,2));
theta   = 180/pi*atan2(-Rtot(1,3),Rtot(3,3));
chi     = asind(-Rtot(2,3));

%% print diagnostics and output

% fprintf(['\n\n\nScattering vector in lab coordinates:',... 
% '\nG(h,k,l) = [%f, %f, %f]\n'], Rtot*Ghkl(:))

fhkl2 = abs(Fhkl(geometry, hkl))^2;

fprintf(['HKL\tThetaB\t|Phi\tTheta\tChi\tMu\tDelta|\t Fhkl^2 \n'])
fprintf(['%2.1f %2.1f %2.1f\t%2.3f\t%f\t'], hkl, thetaB)
fprintf(['%2.3f  %2.3f  %2.3f  %2.3f  %2.3f  %2.1f\n\n'], phi, theta, chi, mu, del, fhkl2)
% n1=Rsample*[0,1,0]';
% alph1=acosd(k0'*Rtot*n1/norm(k0))-90; 

end