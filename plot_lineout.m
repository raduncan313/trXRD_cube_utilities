function [f, d] = plot_lineout(d)
    imsum = sum(d.off.imgs, 3);
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
    
    lineout.sig = (d.on.imgs(y, x, :)./mean(d.off.imgs(y, x, :), 3)); %fix with sub2ind
    lineout.pixels = [y x];

    f = figure
    ax = subplot(1,3,1);
    imagesc(imsum);
    hold on
    set(gca, 'colorscale', 'log')
    plot(x, y, 'Color', 'red', 'LineWidth', 2)

    ax = subplot(1,3,2);
    
    lineout.sig = (d.on.imgs(y, x, :)./mean(d.off.imgs(y, x, :), 3));
    lineout.pixels = [y x];
    d.lineout = lineout;
end