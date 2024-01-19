function [fs, d] = plot_rois(d, n)

    d = clear_masks(d);
    f2 = figure;
    hold on

    for ii = 1:n
        if ii == 1
            [d, f1] = roi_mask(d);
        else
            [d, f1] = roi_mask(d, f1);
        end

        figure(f2);
        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 1.5);
    end

    xlabel('Time (ps)')
    ylabel('Intensity (norm.)')
    fs = {f1, f2};
end