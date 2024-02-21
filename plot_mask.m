function f = plot_mask(d, num)
    f = figure
    subplot(1,2,1)
    imagesc(d.masks{num}.mask)
    title('Mask')
    subplot(1,2,2)
    imagesc(sum(d.off.imgs, 3))
    title('Summed off-image')
end