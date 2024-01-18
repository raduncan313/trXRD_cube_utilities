function [data1d, qs, cx, cy, xi, yi] = qslineout(geometry,data,varargin)

nqy = size(data,2);
nqz = size(data,1);
geometry.imageNy = nqy;
geometry.imageNz = nqz;
% bg=geometry.primvects';

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;

[Q, ~] = plotdata(geometry, data, varargin{1});
[cx, cy, c, xi, yi] = improfile;

Q2d = reshape(reshape(Q,[nqy*nqz, 3])*Rot, [nqz, nqy, 3]);
data1d    = interp2(data, cx, cy);
Q1d_x = interp2(Q2d(:,:,1), cx, cy);
Q1d_y = interp2(Q2d(:,:,2), cx, cy);
Q1d_z = interp2(Q2d(:,:,3), cx, cy);
qs = [Q1d_x Q1d_y Q1d_z];
1;
end
