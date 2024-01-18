function write_sigs_to_csv(ds, folder)
    for ii = 1:length(ds)
        d = ds{ii};
        write_sig_to_csv(d, folder);
    end
end

function write_sig_to_csv(d, folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    savedir = sprintf('%1$s%2$s.csv', folder, d.info);
    savedir = strrep(savedir, ' ', '_');
    writematrix([d.on.scan_var, d.sig], savedir);
end