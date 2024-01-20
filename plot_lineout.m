function [f, d] = plot_lineout(d)
    imsum = sum(d.off.imgs, 3);
    f = figure;
    imagesc(imsum);
    hold on
    set(gca, 'colorscale', 'log');
    title('Select lineout')
    [x, y, ~] = improfile;
    x = round(x);
    y = round(y);
    plot(x, y, 'Color', 'red', 'LineWidth', 2)

    lineout.sig = (d.on.imgs(x, y, :)./mean(d.off.imgs(x,y,:), 3));
    lineout.pixels = [x y];
    %Need to finish
end