function hkl_scatter = anglescan_hklmap(cube, ROI_lims, scan_angle,...
    scan_range, threshs, geometry)

% function hkl_scatter = anglescan_hklmap(cube, ROI_lims, scan_angle,...
%     scan_range, threshs, geometry)
%
%     This function takes a (huber) angle scan cube and returns an
%     hkl-mapping for each pixel
% 
%     INPUTS: cube: A 3D array containing the angle scan data. First index
%     is horizontal, second index is vertical, third index is 'scan_range'
%     index. Care must be taken regarding the orientation of the detector
%     with respect to the indexing of 'cube'. The convention of this
%     function is the standard xdscode convention where the x-axis is the
%     incident x-ray direction, the z-axis points vertically up, and the
%     y-axis completes the right-handed coordinate system. When the
%     detector angles 'delta' and 'nu' are zero the incident x-rays point
%     directly at the detector. If I am at the sample facing the detector
%     in this configuration, the relationship between the detector
%     orientation and the cube indices is given by the following:
%         - values 'cube(1,1,:)' are at the top-left - values
%         'cube(end,1,:)' are at the top-right - values 'cube(1,end,:)' are
%         at the bottom-left - values 'cube(end,end,:)' are at the
%         bottom-right
%     That is, in this configuration increasing the horizontal index takes
%     you along -y, and increasing the vertical index takes you along -z.
% 
%     ROI_lims: a two-element cell with each element itself containing a
%     two-element array defining the boundaries of the ROI. The first index
%     of the cell is horizontal, the second index is vertical. The first
%     index of each of the arrays is the lower value, and the second is the
%     higher value.
% 
%     scan_angle: a char that denotes the angle being scanned. Valid values
%     are 'theta', 'phi', 'chi', and 'mu'.
% 
%     scan_range: a 1D array with the *deviational* angle values of the
%     scan about the given value of the scan angle in the geometry
%     structure.
% 
%     threshs: 2-element array with low (first index) and high (second
%     index) intensity thresholding values. Pixels with intensity values
%     outside this range are skipped in generating 'hkl_scatter'.
% 
%     geometry: standard xdscode geometry structure.
% 
%     The detector angles nu and delta are right-handed and left-handed
%     respectively. The handedness of the sample angles can vary depending
%     on 'geometry.rot_matrix', but this function assumes the general form
%     of 'huber_sixcircle_matrix' or 'huber_fourcircle_matrix'.
% 
%     OUTPUT: hkl_scatter: scattered data corresponding to each pixel in
%     scan mapped to a position in reciprocal space (in r.l.u.). Each row
%     is a given pixel in a given orientation. Columns 1, 2, and 3 are h-,
%     k-, and l-values respectively, column 4 is pixel intensity, and
%     column 5 is the corresponding value of `scan_var`. The intensities
%     are scaled by 1/cos(gamma), where gamma is the angle of the scattered
%     x-ray wavevector relative to the detector surface normal, to yield
%     values proportional to the x-ray fluence incident on the pixel.

    tic
    % Detector distance and pixel size
    L_det = geometry.detector.det_dist; % in millimeters
    l_pix_h = geometry.detector.det_size_horz / geometry.detector.det_pixels_horz; % in millimeters
    l_pix_v = geometry.detector.det_size_vert / geometry.detector.det_pixels_vert; % in millimeters
    l_pix = [l_pix_h, l_pix_v]';

    ROI_lim_h = ROI_lims{1};
    ROI_lim_v = ROI_lims{2};
    ROI_size = (ROI_lim_h(2) - ROI_lim_h(1) + 1)*(ROI_lim_v(2) - ROI_lim_v(1) + 1);
    
    h_ind_c = geometry.beam_center(1);
    v_ind_c = geometry.beam_center(2);
    hvc = [h_ind_c, v_ind_c]';
    
    lambda = geometry.lambda0;
    xunit = [1 0 0]';
    yunit = [0 1 0]';
    zunit = [0 0 1]';
    s_in = (1/lambda)*xunit;
    Rot_S = geometry.rot_matrix;
    Rot_D = @(d,n) rotationmat3D(n,zunit)*rotationmat3D(d,-yunit);
    SamRot = geometry.SamRot;

    if nargin(geometry.rot_matrix) == 4
        fprintf('Assuming sixcircle matrix...\n')
        mattype = 1;
    elseif nargin(geometry.rot_matrix) == 3
        fprintf('Assuming fourcircle matrix...\n')
        mattype = 2;
    else
        error('Invalid `geometry.rot_matrix`: this function only works with sixcircle (4 args) or fourcircle (3 args) matrices.')
    end
    nu = geometry.nu;
    delta = geometry.delta;
    
    num_calc = ROI_size * length(cube(1,1,:));
    hkl_scatter_0 = zeros(num_calc,5);
    p = 1;
    
    for ii = 1:length(cube(1,1,:))
        if strcmp(scan_angle, 'theta')
            theta = geometry.theta + scan_range(ii);
            phi = geometry.phi;
            chi = geometry.chi;
            if mattype == 1
                mu = geometry.mu;
            end
        elseif strcmp(scan_angle, 'phi')
            theta = geometry.theta;
            phi = geometry.phi + scan_range(ii);
            chi = geometry.chi;
            if mattype == 1
                mu = geometry.mu;
            end
        elseif strcmp(scan_angle, 'chi')
            theta = geometry.theta;
            phi = geometry.phi;
            chi = geometry.chi + scan_range(ii);
            if mattype == 1
                mu = geometry.mu;
            end
        elseif strcmp(scan_angle, 'mu')
            if mattype == 2
                error('No `mu` angle for fourcircle matrix.')
            end
            theta = geometry.theta;
            phi = geometry.phi;
            chi = geometry.chi;
            mu = geometry.mu + scan_range(ii);
        end
        
        HV = zeros(2, ROI_size);
        intens = zeros(ROI_size, 1);
        kk = 1;
        for hh = ROI_lim_h(1):ROI_lim_h(2)
            for vv = ROI_lim_v(1):ROI_lim_v(2)
                if (cube(hh,vv,ii) > threshs(1) && cube(hh,vv,ii) < threshs(2))
                    HV(:, kk) = [hh vv]';
                    intens(kk) = cube(hh,vv,ii);
                    kk = kk + 1;  
                end
            end
        end
        npx = kk - 1;
        HV = HV(:,1:npx);
        if isempty(HV)
            continue
        end
        intens = intens(1:npx);        
        R_det = L_det*Rot_D(delta, nu)*xunit;
        R_pix = Rot_D(delta, nu)*(L_det*xunit + [-yunit, -zunit]*((HV - hvc).*l_pix));
        L_pix = vecnorm(R_pix - R_det, 2, 1)';
        gamma = atand(L_pix./L_det);
        ax = cross(repmat(R_det, [1,npx]), R_pix, 1);
        ax = ax./vecnorm(ax, 2, 1);
        
        Rs = zeros(3,3,npx);
        for jj = 1:npx
            gam = gamma(jj);
            axx = ax(:, jj);
            R = rotationmat3D(gam, axx);
            Rs(:,:,jj) = R;
        end
        
        s_out = pagemtimes(Rs, repmat((Rot_D(delta,nu))*s_in, [1,1,npx]));
        s_diff = s_out - s_in;
        if mattype == 1
            s_diff_crystal = pagemldivide(repmat(Rot_S(phi,theta,chi,mu)*SamRot, [1,1,npx]), s_diff);
        else
            s_diff_crystal = pagemldivide(repmat(Rot_S(phi,theta,chi)*SamRot, [1,1,npx]), s_diff);
        end
        hkl = squeeze(permute(pagemtimes(geometry.realvecs, s_diff_crystal), [3 1 2]));        
        hkl_scatter_0(p:p+npx-1,1:3) = hkl;
        hkl_scatter_0(p:p+npx-1,4) = intens./cosd(gamma);
        hkl_scatter_0(p:p+npx-1,5) = scan_range(ii);
        p = p + npx;
        
    end
    hkl_scatter = hkl_scatter_0(1:p-1,:);
    toc
end