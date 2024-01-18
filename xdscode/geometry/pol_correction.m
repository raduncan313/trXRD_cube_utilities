function angcor = pol_correction(geometry)%,image,varargin)
% POL_CORRECTION: Calculate the geometric factor that corrects the
% intensity due to absorption and an additional geometry factor see [M.
% Holt paper and references for details]
% 
% angle_correction = POL_CORRECTION(geometry) returns a matrix of the same
% dimensions as the image size indicated in the input structure geometry.
% This function uses the same logic as DET_KSPACE_PROJ to construct the
% scattering angle and phi.

	det_size_horz = geometry.detector.det_size_horz;
	det_size_vert = geometry.detector.det_size_vert;

    det_pixels_horz = geometry.detector.det_pixels_horz;
	det_pixels_vert = geometry.detector.det_pixels_vert;
    
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
    
    [Ymesh Zmesh] = meshgrid(Y,Z);
    Deno=sqrt(Zmesh.^2 + Ymesh.^2 + D.^2);
    
	Q = zeros([image_Nz, image_Ny, 3]);
    
    Q(:,:,1) = k_0*(D./Deno - 1);           % X-ray incident along +[1 0 0]
    Q(:,:,2) = k_0*( - ( Ymesh./Deno ) );   % RIGHT of image is Y < 0 in the lab frame
    Q(:,:,3) = k_0*(- Zmesh./Deno);         % positive Z is up, but detector local "Z" increases downwards

    ksx = D./Deno;
    twotheta = acos(ksx);
    phi = atan(Ymesh./Zmesh);
    
    angcor = (sin(phi).^2+cos(phi).^2.*cos(twotheta).^2).*cos(twotheta).^3;

    1;
    
