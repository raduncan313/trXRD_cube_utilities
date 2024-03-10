classdef Cube < handle
    properties
        on = struct();
        off = struct();
        lineout = struct();
        runs = {};
        info = {};
        masks = {};
        normed_by_i0 = false;
    end

    methods
        function obj = Cube(cubedir, runs, info)
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
            subplot(1,2,1)
            imagesc(obj.masks{num}.mask)
            title('Mask')
            subplot(1,2,2)
            imagesc(sum(obj.off.imgs, 3))
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