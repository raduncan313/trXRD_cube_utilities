function [f, d] = thresh_and_plot(d, th)
    d = clear_masks(d);
    d = thresh_mask(d, th);

    f = figure
    plot(d.on.scan_var, d.masks{1}.sig, 'LineWidth', 1.5);
    xlabel('Time (ps)')
    ylabel('Intensity (norm.)')
    title(sprintf('%1$s: th = %2$f', d.info, d.masks{1}.th))
end

function d = thresh_mask(d, th)
    imoff = sum(d.off.imgs, 3);
    imoff = imoff/max(imoff(:));
    mask2D = (imoff > th);
    mask3D = repmat(mask2D, [1,1,length(d.off.scan_var)]);
    imgs_mask_on = d.on.imgs.*mask3D;
    imgs_mask_off = d.off.imgs.*mask3D;
    sig = squeeze(sum(imgs_mask_on, [1,2]))./mean(squeeze(sum(imgs_mask_off, [1,2])));
    
    if ~isfield(d, 'masks')
        d.masks = {};
    end
    
    mask.mask = mask2D;
    mask.type = 'thresh';
    mask.th = th;
    mask.sig = sig;
    d.masks{end+1} = mask;
end