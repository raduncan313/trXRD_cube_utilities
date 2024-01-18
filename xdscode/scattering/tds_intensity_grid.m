function [intens, ww, vv] = tds_intensity_grid(geometry, ifc)

ry_to_thz=3289.842;
amu_to_ry=911.4442;

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot     = geometry.rot_matrix(phi0, theta, chi)*SamRot;
% Rot = huber_matrix(phi0, theta, chi)*SamRot;

Q=det_kspace_sphere(geometry);
Qflat=reshape(Q,[],3)*Rot;
nq=size(Qflat,1);
masses = geometry.masses*amu_to_ry;

DD2 = dynmat2(Qflat, ifc, masses);
[ww, vv] = phonon_freqs_vecs(DD2);
intens=tdsI1(geometry, Qflat, ww, vv);
ww = ww*ry_to_thz;
end
