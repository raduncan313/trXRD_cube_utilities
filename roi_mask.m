function d = roi_mask(d)
    imsum = sum(d.off.imgs, 3);
    f = figure;
    title('Select ROI')
    im = imagesc(imsum);
    set(gca, 'colorscale', 'log');
    [mask2D, x, y] = roipoly;
    hold on
    plot(x, y, 'color', 'red')

    mask3D = repmat(mask2D, [1,1,length(d.off.scan_var)]);
    imgs_mask_on = d.on.imgs.*mask3D;
    imgs_mask_off = d.off.imgs.*mask3D;
    sig = squeeze(sum(imgs_mask_on, [1,2]))./mean(squeeze(sum(imgs_mask_off, [1,2])));

    if ~isfield(d, 'masks')
        d.masks = {};
    end
    
    mask.mask = mask2D;
    mask.type = 'ROI';
    mask.roi = [x y];
    mask.sig = sig;
    d.masks{end+1} = mask;
    
    close(f);
end