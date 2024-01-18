function [Q1d, QQ1d] = get_wavevectors_from_path_coord(geometry, img_size, cx, cy)

Ny = img_size(1);
Nx = img_size(2);

geometry.imageNy = Ny;
geometry.imageNz = Nx;

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;

[Q, QQ, ~, ~] = generate_reduced_wavevectors(geometry);

Q1d = vector_interpolate_in_image(Q*Rot,cx,cy,Nx,Ny); % this is in the lab frame, rotated
QQ1d = vector_interpolate_in_image(QQ,cx,cy,Nx,Ny);  % this one is not rotated

end

function Q1d=vector_interpolate_in_image(Q,cx,cy,Nx,Ny)
    Q2d = reshape(Q, [Nx, Ny, 3]);
    Q1d_x = interp2(Q2d(:,:,1), cx, cy);
    Q1d_y = interp2(Q2d(:,:,2), cx, cy);
    Q1d_z = interp2(Q2d(:,:,3), cx, cy);
    Q1d = [Q1d_x Q1d_y Q1d_z];

end
