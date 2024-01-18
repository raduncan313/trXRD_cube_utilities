function Q=det_kspace_sphere(geometry)
% function Q=det_kspace_sphere(geometry)
%
% generates a 2D slice of the Ewald sphere in spherical coordinates

%% detector related parameters and local (2D) coordinate system

	det_size_horz = geometry.detector.det_size_horz;
	det_size_vert = geometry.detector.det_size_vert;

    det_pixels_horz = double(geometry.detector.det_pixels_horz);
	det_pixels_vert = double(geometry.detector.det_pixels_vert);
    
    pix_size_horz = det_size_horz/det_pixels_horz ;
    pix_size_vert = det_size_vert/det_pixels_vert ;
    
    image_Ny = geometry.imageNy;
    image_Nz = geometry.imageNz;
    
    beam_center = geometry.beam_center;
    
    % Local detector coordinate system define the local detector
    % coordinates in image pixels. Origin is at top left pixel
	Y = linspace(0, det_size_horz, image_Ny);	
    Z = linspace(0, det_size_vert, image_Nz);
    
%% stuff in the XYZ (3D) scattering coordinate system

    D = geometry.detector.det_dist;
	k_0 = 1/geometry.lambda0;	%[Angstroms^-1]
    
    % BEAM CENTER:
    % make sure Z and Y are zero at the beam intersection (i.e. [1 0 0])
    Y = Y - beam_center(1)*pix_size_horz;
    Z = Z - beam_center(2)*pix_size_vert;
    
%%
    [Ymesh, Zmesh] = meshgrid(Y,Z);
    Deno=sqrt(Zmesh.^2 + Ymesh.^2 + D.^2);
    
	Q = zeros([image_Nz, image_Ny, 3]);
    
%     [D./Deno, Ymesh./Deno, Zmesh./Deno];

%%
% detRot=rotationmat3D(90, [0,-1,0])
detRot = geometry.detector.detRot;
y = Ymesh(:);
z = Zmesh(:);
x = D*ones(size(z));
rrot = [x,-y,-z]*detRot';

%% plotting the rotated and unrotated pixel grid (sub-sampled)
% clf(); 
% hold all; 
% axis vis3d
% axis equal
% idxs=1:411:size(x,1);
% aa=[x,-y,-z];
% scatter3(rrot(idxs,1),rrot(idxs,2),rrot(idxs,3),'r.')
% scatter3(aa(idxs,1),aa(idxs,2),aa(idxs,3),'b.')
% view(10, 26)

rrot = reshape(rrot, [image_Nz, image_Ny, 3]);

%%
    Q(:,:,1) = k_0*(rrot(:,:,1)./Deno - 1);           % X-ray incident along +[1 0 0]
    Q(:,:,2) = k_0*(( rrot(:,:,2)./Deno ) );   % RIGHT of image is Y < 0 in the lab frame
    Q(:,:,3) = k_0*(rrot(:,:,3)./Deno);         % positive Z is up, but detector local "Z" increases downwards
  
    1;