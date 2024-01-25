function [f, d, inds] = plot_lineout(d)
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

    sz = size(d.on.imgs);
    inds = sub2ind([sz(1), sz(2)], y, x)
    
    imgs_on_2D = reshape(d.on.imgs, [sz(1)*sz(2), sz(3)]);
    imgs_off_2D = reshape(d.off.imgs, [sz(1)*sz(2), sz(3)]);
    sig = imgs_on_2D(inds, :)./mean(imgs_off_2D(inds, :), 2);

    f = figure
    ax = subplot(1,3,1);
    imagesc(imsum);
    hold on
    set(gca, 'colorscale', 'log')
    plot(x, y, 'Color', 'red', 'LineWidth', 2)

    ax = subplot(1,3,2);
    imagesc([d.on.scan_var(1), d.on.scan_var(end)], [], sig);
    xlabel('scan_var')
    
    dt = d.on.scan_var(2) - d.on.scan_var(1);
    freq = linspace(0, 1/dt, length(d.on.scan_var(d.on.scan_var > 0)));
    sig_fft = abs(fft(sig(:,(d.on.scan_var > 0)), [], 2))';
    ax = subplot(1,3,3, 'ydir', 'normal');
    imagesc([], [freq(1), freq(end)], sig_fft)
    set(ax, 'ydir', 'normal');

    lineout.sig = sig;
    lineout.pixels = [y, x];
    lineout.fft = sig_fft;
    d.lineout = lineout;
end