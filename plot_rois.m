function [f, d] = plot_rois(d, n)

    d = clear_masks(d);
    for ii = 1:n
        if ii == 1
            [d, f1] = roi_mask(d);
        else
            [d, f1] = roi_mask(d, f1);
        end

        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 1.5);
    end
    
    close(f1)
    f = figure;
    subplot(1,2,1);
    imsum = sum(d.off.imgs, 3);
    imagesc(imsum);
    set(gca, 'colorscale', 'log');
    hold on
    
    for ii = 1:n
        m = d.masks{ii};
        plot(m.roi(:,1), m.roi(:,2), 'linewidth', 2);
    end
    
    subplot(1,2,2);
    hold on
    
    for ii = 1:n
        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 1.5);
    end
    
    xlabel('Time (ps)')
    ylabel('Intensity (norm.)')
end

function [d, f] = roi_mask(d, varargin)

    if ~isfield(d, 'masks')
        d.masks = {};
    end
    
    imsum = sum(d.off.imgs, 3);
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

    mask3D = repmat(mask2D, [1,1,length(d.off.scan_var)]);
    imgs_mask_on = d.on.imgs.*mask3D;
    imgs_mask_off = d.off.imgs.*mask3D;
    sig = squeeze(sum(imgs_mask_on, [1,2]))./mean(squeeze(sum(imgs_mask_off, [1,2])));
    
    mask.mask = mask2D;
    mask.type = 'ROI';
    mask.roi = [x y];
    mask.sig = sig;
    d.masks{end+1} = mask;
end