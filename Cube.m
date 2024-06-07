classdef Cube < matlab.mixin.Copyable
    properties
        on = struct();
        off = struct();
        lineout = struct();
        runs = {};
        info = {};
        masks = {};
        normed_by_i0 = false;
        geometry = [];
        hkl_scatter = [];
        scan_var_name = '';
    end

    methods
        function obj = Cube(cubedir, runs, scan_var_name, info)
            function d = read_cube(cubedir, run)
                d.runs = {run};
                fname_on = sprintf('%1$s/run%2$04d_on.h5', cubedir, run);
                d.on.scan_var = h5read(fname_on, '/scan_var');
                d.on.imgs = h5read(fname_on, '/imgs');
                d.on.i0 = h5read(fname_on, '/i0');
                d.on.bin_counts = h5read(fname_on, '/bin_counts');    
            
                fname_off = sprintf('%1$s/run%2$04d_off.h5', cubedir, run);
                d.off.scan_var = h5read(fname_off, '/scan_var');
                d.off.imgs = h5read(fname_off, '/imgs');
                d.off.i0 = h5read(fname_off, '/i0');
                d.off.bin_counts = h5read(fname_off, '/bin_counts');    
            end

            for ii = 1:length(runs)
                run = runs{ii};
                d = read_cube(cubedir, run);
                if ii == 1
                    on = d.on;
                    off = d.off;
                else
                    on.imgs = on.imgs + d.on.imgs;
                    on.bin_counts = on.bin_counts + d.on.bin_counts;
                    on.i0 = on.i0 + d.on.i0;
                    
                    off.imgs = off.imgs + d.off.imgs;
                    off.bin_counts = off.bin_counts + d.off.bin_counts;
                    off.i0 = off.i0 + d.off.i0;
                end
            end

            obj.on = on;
            obj.off = off;
            obj.runs = runs;
            obj.info = info;
            obj.scan_var_name = scan_var_name;
        end

        function add_geometry(obj, geometry)
            obj.geometry = geometry;
        end

        function transpose(obj)
            obj.clear_masks();
            obj.clear_lineout();
            obj.clear_geometry();
            obj.on.imgs = permute(obj.on.imgs, [2,1,3]);
            obj.off.imgs = permute(obj.off.imgs, [2,1,3]);
        end

        function flip(obj, dim)
            obj.clear_masks();
            obj.clear_lineout();
            obj.clear_geometry();
            obj.on.imgs = flip(obj.on.imgs, dim);
            obj.off.imgs = flip(obj.off.imgs, dim);
        end

        function norm_i0(obj)
            if ~obj.normed_by_i0
                sz = size(obj.on.imgs);
                obj.on.imgs = obj.on.imgs./permute(repmat(obj.on.i0, [1, sz(1), sz(2)]), [2,3,1]);
                obj.off.imgs = obj.off.imgs./permute(repmat(obj.off.i0, [1, sz(1), sz(2)]), [2,3,1]);
                obj.normed_by_i0 = true;
            else
                fprintf('Already normalized by i0, did not re-normalize.\n')
            end
        end

        function subtract_t0(obj, t0)
            obj.on.scan_var = obj.on.scan_var - t0;
            obj.off.scan_var = obj.off.scan_var - t0;
        end

        function clear_masks(obj)
            obj.masks = {};
        end

        function clear_lineout(obj)
            obj.lineout = struct();
        end

        function clear_geometry(obj)
            obj.geometry = [];
        end

        function f = thresh_and_plot(obj, th)
            obj.clear_masks();
            obj.thresh_mask(th);
        
            f = figure
            plot(obj.on.scan_var, obj.masks{1}.sig, 'LineWidth', 1.5);
            xlabel('Time (ps)')
            ylabel('Intensity (norm.)')
            title(sprintf('%1$s: th = %2$f', obj.info, obj.masks{1}.th))
        end

        function rebin(obj, binint)
            obj.clear_masks();
            obj.clear_lineout();

            scan_var_on_0 = obj.on.scan_var;
            imgs_on_0 = obj.on.imgs;
            i0_on_0 = obj.on.i0;
            bin_counts_on_0 = obj.on.bin_counts;
            
            scan_var_off_0 = obj.off.scan_var;
            imgs_off_0 = obj.off.imgs;
            i0_off_0 = obj.off.i0;
            bin_counts_off_0 = obj.off.bin_counts;
        
            len1 = length(imgs_on_0(:,1,1));
            len2 = length(imgs_on_0(1,:,1));
            len3_0 = length(imgs_on_0(1,1,:));
            len3 = floor(len3_0/binint);
        
            scan_var_on = zeros(len3,1);    
            imgs_on = zeros(len1,len2,len3);
            i0_on = zeros(len3,1);
            bin_counts_on = zeros(len3,1);
            
            scan_var_off = zeros(len3,1);
            imgs_off = zeros(len1,len2,len3);
            i0_off = zeros(len3,1);
            bin_counts_off = zeros(len3,1);
        
            hh = 1;
            for ii = 1:len3
                imb_on = zeros(len1,len2);
                i0b_on = 0;
                bcb_on = 0;
                
                imb_off = zeros(len1,len2);
                i0b_off = 0;
                bcb_off = 0;
                
                for jj = 1:binint
                    imb_on = imb_on + imgs_on_0(:,:,hh);
                    i0b_on = i0b_on + i0_on_0(hh);
                    bcb_on = bcb_on + bin_counts_on_0(hh);
                    
                    imb_off = imb_off + imgs_off_0(:,:,hh);
                    i0b_off = i0b_off + i0_off_0(hh);
                    bcb_off = bcb_off + bin_counts_off_0(hh);
                    
                    hh = hh + 1;
                end
                
                imgs_on(:,:,ii) = imb_on;
                i0_on(ii) = i0b_on;
                bin_counts_on(ii) = bcb_on;
                
                imgs_off(:,:,ii) = imb_off;
                i0_off(ii) = i0b_off;
                bin_counts_off(ii) = bcb_off;
                
                scan_var_on(ii) = mean([scan_var_on_0(hh-1), scan_var_on_0(hh-binint)]);
                scan_var_off(ii) = mean([scan_var_off_0(hh-1), scan_var_off_0(hh-binint)]);
            end
            
            obj.on.scan_var = scan_var_on;
            obj.on.imgs = imgs_on;
            obj.on.i0 = i0_on;
            obj.on.bin_counts = bin_counts_on;
            
            obj.off.scan_var = scan_var_off;
            obj.off.imgs = imgs_off;
            obj.off.i0 = i0_off;
            obj.off.bin_counts = bin_counts_off;
        end

        function f = plot_rois(obj, n)        
            obj.clear_masks();
            for ii = 1:n
                if ii == 1
                    f1 = obj.roi_mask();
                else
                    f1 = obj.roi_mask(f1);
                end        
                plot(obj.on.scan_var, obj.masks{ii}.sig, 'linewidth', 1.5);
            end            
            close(f1)
            f = figure;
            subplot(1,2,1);
            imsum = sum(obj.off.imgs, 3);
            imagesc(imsum);
            set(gca, 'colorscale', 'log');
            hold on            
            for ii = 1:n
                m = obj.masks{ii};
                plot(m.roi(:,1), m.roi(:,2), 'linewidth', 2);
            end            
            subplot(1,2,2);
            hold on            
            for ii = 1:n
                plot(obj.on.scan_var, obj.masks{ii}.sig, 'linewidth', 1.5);
            end            
            xlabel('scan\_var')
            ylabel('Intensity (norm.)')
        end

        function write_sigs_to_csv(obj, folder)
            if ~exist(folder, 'dir')
                mkdir(folder);
            end            
            for ii = 1:length(obj.masks)
                savedir = sprintf('%1$s/%2$s_mask%3$d.csv', folder, obj.info, ii);
                savedir = strrep(savedir, ' ', '_');
                writematrix([obj.on.scan_var, obj.masks{ii}.sig], savedir);
            end
        end

        function plot_sigs_to_axis(obj, ax, leg)        
            axes(ax);
            hold on
            leg_str = leg.String;            
            for ii = 1:length(obj.masks)
                plot(obj.on.scan_var, obj.masks{ii}.sig, 'linewidth', 1.5)
                leg_str{end+1} = sprintf('%1$s mask %2$d', obj.info, ii);
            end            
            set(leg, 'string', leg_str)
        end

        function f = plot_mask(obj, num)
            f = figure
            ax = subplot(1,2,1)
            imagesc(obj.masks{num}.mask);
            title('Mask')
            ax = subplot(1,2,2);
            imagesc(sum(obj.off.imgs, 3));
            set(ax, 'colorscale', 'log');
            title(sprintf('Summed off-image: %1$s', obj.info))
        end
        
        function f = plot_lineout(obj)
            imsum = sum(obj.off.imgs, 3);
            f1 = figure;
            imagesc(imsum);
            hold on
            set(gca, 'colorscale', 'log');
            title('Select lineout')
            [x, y, ~] = improfile;
            x = round(x);
            y = round(y);
            plot(x, y, 'Color', 'red', 'LineWidth', 2)
            close(f1)
        
            sz = size(obj.on.imgs);
            inds = sub2ind([sz(1), sz(2)], y, x)
            
            imgs_on_2D = reshape(obj.on.imgs, [sz(1)*sz(2), sz(3)]);
            imgs_off_2D = reshape(obj.off.imgs, [sz(1)*sz(2), sz(3)]);
            sig = imgs_on_2D(inds, :)./mean(imgs_off_2D(inds, :), 2);
        
            f = figure
            ax = subplot(1,3,1);
            imagesc(imsum);
            hold on
            set(gca, 'colorscale', 'log')
            plot(x, y, 'Color', 'red', 'LineWidth', 2)
        
            ax = subplot(1,3,2);
            imagesc([obj.on.scan_var(1), obj.on.scan_var(end)], [], sig);
            xlabel('scan\_var')
            
            dt = obj.on.scan_var(2) - obj.on.scan_var(1);
            freq = linspace(0, 1/dt, length(obj.on.scan_var(obj.on.scan_var > 0)));
            sig_fft = abs(fft(sig(:,(obj.on.scan_var > 0)), [], 2))';
            ax = subplot(1,3,3, 'ydir', 'normal');
            imagesc([], [freq(1), freq(end)], sig_fft)
            set(ax, 'ydir', 'normal');
        
            lineout.sig = sig;
            lineout.pixels = [y, x];
            lineout.fft = sig_fft;
            lineout.inds = inds
            obj.lineout = lineout;
        end

        function f = plotdatacube(obj, type, caxis, varargin)
            if strcmp(type, 'on')
                data = obj.on.imgs;
            elseif strcmp(type, 'off')
                data = obj.off.imgs;
            elseif strcmp(type, 'ratio')
                data = obj.on.imgs./mean(obj.off.imgs, 3);
            else
                error('Invalid `type` parameter -- must be `on`, `off`, or `ratio`.')
            end

            if nargin == 4
                if strcmp(varargin{1}, 'log')
                    data = log10(data);
                end
            end
            plotdatacube(obj.on.scan_var, data, obj.geometry, caxis, []);
            f = gcf;
        end

        function f = hklmap(obj, ROI_lims, threshs)
            if isempty(obj.geometry)
                error('Need a `geometry` struct for hkl mapping.')
            end
            
            L_det = obj.geometry.detector.det_dist;
            l_pix_h = obj.geometry.detector.det_size_horz / obj.geometry.detector.det_pixels_horz;
            l_pix_v = obj.geometry.detector.det_size_vert / obj.geometry.detector.det_pixels_vert;
            l_pix = [l_pix_h, l_pix_v]';
        
            ROI_lim_h = ROI_lims{1};
            ROI_lim_v = ROI_lims{2};
            ROI_size = (ROI_lim_h(2) - ROI_lim_h(1) + 1)*(ROI_lim_v(2) - ROI_lim_v(1) + 1);
            
            h_ind_c = obj.geometry.beam_center(1);
            v_ind_c = obj.geometry.beam_center(2);
            hvc = [h_ind_c, v_ind_c]';
            
            lambda = obj.geometry.lambda0;
            xunit = [1 0 0]';
            yunit = [0 1 0]';
            zunit = [0 0 1]';
            s_in = (1/lambda)*xunit;
            Rot_S = obj.geometry.rot_matrix;
            Rot_D = @(d,n) rotationmat3D(n,zunit)*rotationmat3D(d,-yunit);
            SamRot = obj.geometry.SamRot;
        
            if nargin(obj.geometry.rot_matrix) == 4
                fprintf('Assuming sixcircle matrix...\n')
                mattype = 1;
            elseif nargin(obj.geometry.rot_matrix) == 3
                fprintf('Assuming fourcircle matrix...\n')
                mattype = 2;
            else
                error('Invalid `geometry.rot_matrix`: this function only works with sixcircle (4 args) or fourcircle (3 args) matrices.')
            end
            nu = obj.geometry.nu;
            delta = obj.geometry.delta;
            
            num_calc = ROI_size * length(obj.off.imgs(1,1,:));
            hkl_scatter_0 = zeros(num_calc,5);
            p = 1;
            
            for ii = 1:length(obj.off.imgs(1,1,:))
                if strcmp(obj.scan_var_name, 'theta')
                    theta = obj.off.scan_var(ii);
                    phi = obj.geometry.phi;
                    chi = obj.geometry.chi;
                    if mattype == 1
                        mu = obj.geometry.mu;
                    end
                elseif strcmp(obj.scan_var_name, 'phi')
                    theta = obj.geometry.theta;
                    phi = obj.off.scan_var(ii);
                    chi = obj.geometry.chi;
                    if mattype == 1
                        mu = obj.geometry.mu;
                    end
                elseif strcmp(obj.scan_var_name, 'chi')
                    theta = obj.geometry.theta;
                    phi = obj.geometry.phi;
                    chi = obj.off.scan_var(ii);
                    if mattype == 1
                        mu = obj.geometry.mu;
                    end
                elseif strcmp(obj.scan_var_name, 'mu')
                    if mattype == 2
                        error('No `mu` angle for fourcircle matrix.')
                    end
                    theta = obj.geometry.theta;
                    phi = obj.geometry.phi;
                    chi = obj.geometry.chi;
                    mu = obj.off.scan_var(ii);
                else
                    error('Invalid `scan_var_name`.')
                end
                
                HV = zeros(2, ROI_size);
                intens = zeros(ROI_size, 1);
                kk = 1;
                for hh = ROI_lim_h(1):ROI_lim_h(2)
                    for vv = ROI_lim_v(1):ROI_lim_v(2)
                        if (obj.off.imgs(hh,vv,ii) > threshs(1) && obj.off.imgs(hh,vv,ii) < threshs(2))
                            HV(:, kk) = [hh vv]';
                            intens(kk) = obj.off.imgs(hh,vv,ii);
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
                hkl = squeeze(permute(pagemtimes(obj.geometry.realvecs, s_diff_crystal), [3 1 2]));        
                hkl_scatter_0(p:p+npx-1,1:3) = hkl;
                hkl_scatter_0(p:p+npx-1,4) = intens./cosd(gamma);
                hkl_scatter_0(p:p+npx-1,5) = obj.off.scan_var(ii);
                p = p + npx;
                
            end
            obj.hkl_scatter = hkl_scatter_0(1:p-1,:);
            f = figure
            scatter3(obj.hkl_scatter(:,1),obj.hkl_scatter(:,2),...
                obj.hkl_scatter(:,3),...
                10*ones(length(obj.hkl_scatter(:,1)),1),...
                obj.hkl_scatter(:,5),'filled');
            xlabel('h')
            ylabel('k')
            zlabel('l')
            axis vis3d
        end

        function f = LPSVD_sig(obj, num, L, rat, t_start)
            if ~strcmp(obj.scan_var_name, 'delay')
                error('Invalid run type -- can only use LPSVD with `delay` scans.')
            end
            trunc = obj.on.scan_var >= t_start;
            t = obj.on.scan_var(trunc);
            x = obj.masks{num}.sig(trunc)
            dt = t(2) - t(1);
            N = length(x);
            M = floor(rat*N);
            X = zeros(N - M, M);
            
            for i=1:M
               X(:,i)=x(i+1:i+N-M);
            end
            
            [U, S, V] = svd(X);
            Sinv = 1./S';
            Sinv(isinf(Sinv)) = 0;
            Sinv(isnan(Sinv)) = 0;
            
            for ii = L+1:min(N-M, M)
                Sinv(ii,ii) = 0;
            end
        
            a = [1; -V*Sinv*U'*x(1:N-M)];
            rts = roots(a);
            rts_srt = flip(sort(rts));
            rts_srt_trunc = rts_srt(1:rank(S));
            bs = real(log(rts_srt_trunc))/dt;
            ws = imag(log(rts_srt_trunc))/dt;
            
            bs = bs(ws >= 0);
            ws = ws(ws >= 0);
            
            ws = ws(bs >= 0);
            bs = bs(bs >= 0);
            
            [ws, I] = sort(ws);
            bs = bs(I);
            
            K = length(ws);
            Y = zeros(N, 2*K + 1);
            
            for ii = 1:K
                b = bs(ii);
                w = ws(ii);
                ycos = exp(-b*t).*cos(w*t);
                ysin = exp(-b*t).*sin(w*t);
                Y(:,2*ii-1) = ycos;
                Y(:,2*ii) = ysin;
            end
            
            Y(:,2*K+1) = 1;
            coefs = Y \ x;
            
            phs = zeros(K,1);
            amps = zeros(K,1);
            for ii = 1:K
                c_cos = coefs(2*ii-1);
                c_sin = coefs(2*ii);
                phs(ii) = atan2(c_sin, c_cos);
                amps(ii) = sqrt(c_sin^2 + c_cos^2);
            end
            
            z = zeros(N,1);
            y = zeros(N, K+1);
            for ii = 1:K
                b = bs(ii);
                ph = phs(ii);
                amp = amps(ii);
                w = ws(ii);
                yy = amp*exp(-b*t).*cos(w*t - ph);
                z = z + yy;
                y(:,ii) = yy;
            end
            y(:,end) = coefs(end);
            z = z + coefs(end);
            
            sol.w = ws;
            sol.b = bs;
            sol.ph = phs;
            sol.amp = amps;
            sol.t = t;
            sol.x = x;
            sol.y = y;
            sol.z = z;
            sol.ssd = sum((x - z).^2);
            sol.L = L;
            sol.rat = rat;
            obj.masks{num}.LPSVD = sol;

            f = figure;
            ax = subplot(1,2,1);
            hold on
            plot(obj.on.scan_var, obj.masks{num}.sig, 'linewidth', 2);
            plot(sol.t, sol.z, 'linewidth', 2);

            ax = subplot(1,2,2);
            hold on
            for ii = 1:K+1
                plot(sol.t, sol.y(:,ii), 'LineWidth', 2);
            end
        end

        function f = auto_signal(obj, numcomponents, epsilon, minpts, frac)
            sz = size(obj.on.imgs);
            on_flat = reshape(obj.on.imgs, [], sz(3));
            on_flat_sc = reshape(obj.on.imgs, [], sz(3));
            mean_on = mean(on_flat_sc, 2);
            std_on = std(on_flat_sc, 0, 2);
            on_flat_sc = (on_flat_sc - mean_on)./std_on;
            on_flat_sc(isnan(on_flat_sc)) = 0;
            on_flat_sc(isinf(on_flat_sc)) = 0;

            off_flat = reshape(obj.off.imgs, [], sz(3));

            [~, X_on, ~] = pca(on_flat_sc, 'NumComponents', numcomponents);

            rng(42);
            [Mdl, has_signal] = ocsvm(X_on, 'contaminationfraction', frac,...
                'kernelscale', 'auto');

            inds_1D = find(has_signal == true);
            [x, y] = ind2sub([sz(1), sz(2)], inds_1D);
            X_db = [x y];
            idx = dbscan(X_db, epsilon, minpts);
            idx_unq = unique(idx);
            idx_unq = idx_unq(idx_unq ~= -1);
            sigs_clsts = cell(1,length(idx_unq));
            inds_clsts = cell(1,length(idx_unq));

            for ii = 1:length(idx_unq)
                inds_clst = X_db(idx == idx_unq(ii), :);
                jj1D = sub2ind([sz(1), sz(2)], inds_clst(:,1), inds_clst(:,2));
                sig_clst = mean(on_flat(jj1D,:), 1)./mean(off_flat(jj1D,:), 'all');
                sigs_clsts{ii} = sig_clst;
                inds_clsts{ii} = inds_clst;
            end
            
            f = figure;
            ax = subplot(1,2,1);
            imagesc(sum(obj.on.imgs, 3))
            hold on
            set(ax, 'colorscale', 'log')
            for ii = 1:length(inds_clsts)
                plot(inds_clsts{ii}(:,2), inds_clsts{ii}(:,1), '.', 'MarkerSize', 10);
            end

            ax = subplot(1,2,2);
            hold on
            for ii = 1:length(sigs_clsts)
                plot(obj.on.scan_var, sigs_clsts{ii} + ii, 'linewidth', 1.5);
            end
        end
    end

    methods (Access = private)
        function thresh_mask(obj, th)
            imoff = sum(obj.off.imgs, 3);
            imoff = imoff/max(imoff(:));
            mask2D = (imoff > th);
            mask3D = repmat(mask2D, [1,1,length(obj.off.scan_var)]);
            imgs_mask_on = obj.on.imgs.*mask3D;
            imgs_mask_off = obj.off.imgs.*mask3D;
            sig = squeeze(sum(imgs_mask_on, [1,2]))./mean(squeeze(sum(imgs_mask_off, [1,2])));

            mask.mask = mask2D;
            mask.type = 'thresh';
            mask.th = th;
            mask.sig = sig;
            obj.masks{end+1} = mask;
        end

        function f = roi_mask(obj, varargin)
            imsum = sum(obj.off.imgs, 3);
            if nargin == 1
                f = figure;
                im = imagesc(imsum);
                hold on
                set(gca, 'colorscale', 'log');
                title('Select ROI')
            else
                f = varargin{1};
                figure(f);
                hold on
            end
            
            [mask2D, x, y] = roipoly;
            plot(x, y, 'linewidth', 2);
        
            mask3D = repmat(mask2D, [1,1,length(obj.off.scan_var)]);
            imgs_mask_on = obj.on.imgs.*mask3D;
            imgs_mask_off = obj.off.imgs.*mask3D;
            sig = squeeze(sum(imgs_mask_on, [1,2]))./mean(squeeze(sum(imgs_mask_off, [1,2])));
            
            mask.mask = mask2D;
            mask.type = 'ROI';
            mask.roi = [x y];
            mask.sig = sig;
            obj.masks{end+1} = mask;
        end
    end
end