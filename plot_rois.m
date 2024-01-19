function [f, d] = plot_rois(d, n)

    d = clear_masks(d);

    for ii = 1:n
        if ii == 1
            [d, f] = roi_mask(d);
        else
            [d, f] = roi_mask(d, f);
        end

        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 1.5);
    end
    
    close(f)
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