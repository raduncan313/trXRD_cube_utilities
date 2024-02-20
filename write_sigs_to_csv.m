function write_sigs_to_csv(d, folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    
    for ii = 1:length(d.masks)
        savedir = sprintf('%1$s/%2$s_mask%3$d.csv', folder, d.info, ii);
        savedir = strrep(savedir, ' ', '_');
        writematrix([d.on.scan_var, d.masks{ii}.sig], savedir);
    end
end