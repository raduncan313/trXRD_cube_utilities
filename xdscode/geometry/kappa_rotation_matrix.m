function Rotmat = kappa_rotation_matrix(phi, kappa, omega, mu, kappa_arm_angle)
% returns the rotation matrix of a kappa diffractometer in our coordinate
% system where:
% X is along the x-ray beam
% 
% the axes for all angles set to zero are:

% kappa_arm_angle = 60;

omega_axis = [0 1 0];
% zaxis = [0,0,1];
kappa_axis = rotationmat3D(kappa_arm_angle, [1,0,0])*[0,1,0]';
% 
% kappa_axis = rotationmat3D(kappa_arm_angle, zaxis)*[0,1,0]';
kappamat = rotationmat3D(kappa, kappa_axis);
Rotmat = rotationmat3D(mu, [0,0,1])*rotationmat3D(omega, omega_axis)*kappamat*rotationmat3D(phi, omega_axis);

end

% according to:
% https://7id.xray.aps.anl.gov/calculators/kappa.html
% we should have:
% Theta Chi Phi         Theta  Kappa   Phi
% 0     10  0           -4.21   13.066  -4.21
% 0,    25  0           -10.721 32.824  -10.721
% 

% 15    25  67          25.721 32.824  77.721
