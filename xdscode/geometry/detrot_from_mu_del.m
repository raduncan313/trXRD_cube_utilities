function R=detrot_from_mu_del(mu, del)
% DETROT_FROM_MU_DEL finds the rotation matrix from the `mu` and `del`
% angles for the detector in the huber diffractometer geometry
% -18.273  20.803
minusy=[0,-1,0]; % -y axis in the lab frame
R1 = rotationmat3D(del, minusy);
R2 = rotationmat3D(mu, [0,0,1]);
R = R2*R1;
end
