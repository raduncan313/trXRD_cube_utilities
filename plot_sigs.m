function f = plot_sigs(d)
    masks = d.masks;
    f = figure;
    title(d.info)
    hold on
    leg = {};
    for ii = 1:length(masks)
        plot(d.on.scan_var, d.masks{ii}.sig, 'linewidth', 2);
        leg{end + 1} = sprintf('mask %1$d', ii);
    end
    legend(leg);
end