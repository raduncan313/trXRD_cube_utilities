%%
addpath(genpath('~/Dropbox/xdscode/'))

load sto_cubic_110K.mat

geometry.lambda0=12.398/9;
at = geometry.realvecs';
bg = geometry.primvects';

ry_to_thz=3289.842;
amu_to_ry=911.4442;

exit_angle= 1;
sample_normal = [0,0,1];
hkl = [1,2,3];

k0 = 1/geometry.lambda0*[1;0;0];

%% test inputs to find_rotation_Bragg...
Rot = find_rotation_Bragg_condition(geometry, hkl);
Rot = find_rotation_Bragg_condition(geometry, hkl, 'sampleNormal', [1,0,0],'Angle',1,'grazingExit', true,'initialguess',11);
Rot = find_rotation_Bragg_condition(geometry, hkl, 'sampleNormal', [1,0,0],'Angle',1,'grazingincidence', true,'initialguess',11);
Rot = find_rotation_Bragg_condition(geometry, hkl,'initialguess',-10);
Rot = find_rotation_Bragg_condition(geometry, hkl, 'Angle',18);

%% Find symmetric Bragg and parametrize the rotation in terms of the xpp angles 
hkl = [0,0,3];
sample_normal = [0,0,1];

Rot = find_rotation_Bragg_condition(geometry, hkl,'sampleNormal', sample_normal);
ff = @(angs) huber_matrix(angs(1), angs(2), angs(3));
find_angles_from_rotation(ff, Rot, [0, 0, 0])

% Here the incidence angle should be the Bragg angle:
fprintf('Angle between k0 and sample normal: %2.2f\n', asind(k0'*Rot*sample_normal(:)/norm(k0)))
thB = asind(geometry.lambda0*norm(geometry.primvects'*hkl(:))/2);
fprintf('Should be the same as theta_Bragg: %2.2f\n', thB)

%% Assymetric reflection with fixed *incidence* angle
hkl = [1,0,2];
sample_normal = [0,0,1];
Rot = find_rotation_Bragg_condition(geometry, hkl,'sampleNormal', sample_normal,'grazingIncidence', true,'Angle',1);
ff = @(angs) huber_matrix(angs(1), angs(2), angs(3));
find_angles_from_rotation(ff, Rot, [0, 0, 0])

% the incidence angle is now fixed:
fprintf('Angle between k0 and sample normal: %2.2f\n', asind(k0'*Rot*sample_normal(:)/norm(k0)))

%% Assymetric reflection with fixed *exit* angle
hkl = [1,0,2];
sample_normal = [0,0,1];
Rot = find_rotation_Bragg_condition(geometry, hkl,'sampleNormal', sample_normal,'grazingExit', true,'Angle',1);
ff = @(angs) huber_matrix(angs(1), angs(2), angs(3));
find_angles_from_rotation(ff, Rot, [0, 0, 0])

% the exit angle is:
kp = k0(:) + Rot*bg*hkl(:);
fprintf('Angle between kp and sample normal: %2.2f\n', asind(kp'*Rot*sample_normal(:)/norm(kp)))


%% compare the new find_rotation_Bragg with the older find_rotation_grazing_exit
exit_angle= 1;
sample_normal = [0,1,0];
hkl = [2,1,3];
geometry.rot_matrix=@(a,b,c) huber_sixcircle_matrix(a,b,c,0);
% this is the "old" function we used during SwissFEL beamtime in Dec 2020,
% which returns a rotation matrix `Rot0` 
Rot0=find_rotation_grazing_exit(geometry, sample_normal, hkl, exit_angle, 0)
% we then use `Rot0` to find the angles in the six-circle diffractometer
% geometry:
find_angles_from_rotation(@(angs) huber_sixcircle_matrix(angs(1),angs(2),angs(3),0), Rot0, [100, 20, 60])

% this new function should return an equivalent rotation matrix as
% `find_rotation_grazing_exit` as above:
Rot = find_rotation_Bragg_condition(geometry,hkl,'grazingExit',true,'sampleNormal',sample_normal,'Angle',exit_angle);
% and the angles should be the same as found above:
find_angles_from_rotation(@(angs) huber_sixcircle_matrix(angs(1),angs(2),angs(3),0), Rot, [100, 20, 60])


%% compare with find_angles_outgoink
% note that `find_rotation_Bragg...` finds *a* rotation, not necessarily
% the same rotation as ...outgoink, since any rotation along the incidence
% k-vector satisfied both Bragg and the incidence angle condition. In order
% to compare them we need to allow for an arbitrary rotation along k0
% toghether with the `phi` and `theta` angles in huber_matrix


% find_angles_outgoingk(geometry, [1,2,3], [0,1,0], 10, -33, 0);
% HKL	ThetaB	|Phi	Theta	Chi	Mu	Delta|	 Fhkl^2 
% 1 2 3	41.297	98.794  10.000  -33.000  80.597  37.916  744.1
% clc
sample_normal=[0,0 -1];
hkl = [1,2,3];

the = 20;
geometry.rot_matrix = @huber_matrix;
[phi, theta, chi, ~] = find_angles_outgoingk(geometry,  hkl, sample_normal, the, 0, 0);
Rot0=huber_matrix(phi, theta, chi);
ff = @(angs) huber_matrix(angs(1), angs(2), angs(3));
find_angles_from_rotation(ff, Rot0, [0, 0, 0])

%
Rot = find_rotation_Bragg_condition(geometry, hkl,'sampleNormal', sample_normal,'grazingIncidence',true,'Angle', the, 'initialGuess', 0);
fprintf('Angle between k0 and sample normal: %2.2f\n', asind(k0'*Rot*sample_normal(:)/norm(k0)))

% ff = @(angs) huber_matrix(angs(1), angs(2), angs(3));
ff = @(angs) rotationmat3D(angs(3),[1,0,0])*huber_matrix(angs(1), angs(2), 0);
solff=find_angles_from_rotation(ff, Rot, [0,0,0])

% allowing the extra rotation along 100 is necessary because `Rot` and
% `Rot0` may differ by a rotation along the incident beam.

