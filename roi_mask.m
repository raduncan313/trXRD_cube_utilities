function [d, varargout] = roi_mask(d, varargin)

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
    
    if nargout == 2
        varargout{1} = f;
    else
        close(f);
    end
end