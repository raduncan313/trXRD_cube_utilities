function plot_sigs_to_axis(d, ax, leg)

    axes(ax);
    hold on
    leg_str = leg.String;
    
    for ii = 1:length(d.masks)
        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 1.5)
        leg_str{end+1} = sprintf('%1$s mask %2$d', d.info, ii);
    end
    
    set(leg, 'string', leg_str)
end