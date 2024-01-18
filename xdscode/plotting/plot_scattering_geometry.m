function plot_scattering_geometry(geometry, hkl, sample_normal, sample_inplane)
% Plot the scattering geometry 
% this uses the package `arrow3` from matlabcentral to plot: 
% https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3/

% sample_inplane = [0,1,0];
lambda0=geometry.lambda0;
bg = geometry.primvects';
hkl=hkl(:);
Ghkl = bg*hkl;


k0 = 1/lambda0*[1,0,0]';

% parameters for plotting arrows
arrow_width = 3;
arrow_height= 5;
o = [0,0,0];
X = [1,0,0];
Y = [0,1,0];
Z = [0,0,1];

% get the rotation matrix for the diffraction condition
% Rtot = find_rotation_grazing_exit(geometry, sample_normal, hkl, beta);
% rather use the Huber angles in `geometry` to define the scattering geometry:
Rtot = geometry.rot_matrix(geometry.phi, geometry.theta, geometry.chi)*geometry.SamRot;

%%% the final rotated vectors are:
finalRotGhkl    = Rtot*Ghkl;
finalRotNsample = Rtot*sample_normal(:);
finalkp         = k0 + finalRotGhkl;

hold on

% plot axes
arrow3(o, X, 'r:', arrow_width/3, arrow_height)
arrow3(o, Y, 'g:', arrow_width/3, arrow_height)
arrow3(o, Z, 'b:', arrow_width/3, arrow_height)

% plot `k0`
arrow3(o, k0(:)', 'k2', arrow_width, arrow_height)
% plot `rotGhkl`
arrow3(o, finalRotGhkl(:)', 'm2', arrow_width, arrow_height)
% plot `kp`
arrow3(o, finalkp(:)', 'g2', arrow_width, arrow_height)

% plot the rotated normal `rotNsample`
arrow3(o, 0.5*finalRotNsample(:)', 'r2', arrow_width, arrow_height)

% Make a square starting from `sample_normal` and `sample_inplane` 
v2 = cross(sample_normal(:)', sample_inplane(:)');
vs = [o; sample_inplane(:)'; v2; v2+sample_inplane(:)'];
vs = vs-mean(vs,1) + o; % center at `o`
faces = [1,2,4,3];

p1 = patch('vertices',vs, 'faces',faces);
p1.FaceAlpha=0.8;
p1.FaceColor=[.95,.5,.4];

% find the final rotation axis and angle
[V, e] = eig(Rtot);
e = diag(e);
rot_axis = V(:,imag(e) == 0);
% the trace of Rtot is 1 + 2*cos(th) with th the rotation angle
rot_angle= acosd((trace(Rtot) - 1)/2);

% to fix the sign of the angle (not well-defined in the trace above) do
% this:
if norm(rotationmat3D(rot_angle, rot_axis) - Rtot) > 1e-15
    rot_angle = -rot_angle;
end

rotate(p1, rot_axis, rot_angle, o);

xlabel('x')
ylabel('y')
zlabel('z')

grid
axis tight vis3d 
rotate3d on

view(-90, 7)
drawnow

hold off

end