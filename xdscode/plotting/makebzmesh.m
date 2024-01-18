function bz=makebzmesh(geometry, datasize)
nz = datasize(1);
ny = datasize(2);
% geometry.imageNz = floor(nz/8);
% geometry.imageNy = floor(ny/8);
gnz = geometry.imageNz;
gny = geometry.imageNy;
bg = geometry.primvects';

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;

Q = det_kspace_proj(geometry);
allhkl = findclosesthkl2(bg, reshape(Q,[],3)*inv(Rot)');
allhkl = reshape(allhkl, [gnz, gny, 3]);

% d2=del2(abs(sqrt(2)*allhkl(:,:,1) - pi*allhkl(:,:,2) + exp(1)*allhkl(:,:,3))); 
d2 = imresize(abs(sqrt(2)*allhkl(:,:,1).^2 - pi*allhkl(:,:,2).^2 + exp(1)*allhkl(:,:,3).^2), [nz, ny], 'nearest');
bz = del2(d2);
1;
