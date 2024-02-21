function d_sum = read_cubes(cubedir, runs, info)
    for ii = 1:length(runs)
        run = runs{ii};
        d = read_cube(cubedir, run);
        if ii == 1
            d_sum = d;
        else
            d_sum.on.imgs = d_sum.on.imgs + d.on.imgs;
            d_sum.on.bin_counts = d_sum.on.bin_counts + d.on.bin_counts;
            d_sum.on.i0 = d_sum.on.i0 + d.on.i0;
            
            d_sum.off.imgs = d_sum.off.imgs + d.off.imgs;
            d_sum.off.bin_counts = d_sum.off.bin_counts + d.off.bin_counts;
            d_sum.off.i0 = d_sum.off.i0 + d.off.i0;
        end
    end
    d_sum.runs = runs;
    d_sum.info = info;
end

function d = read_cube(cubedir, run)
    d.runs = {run};
    fname_on = sprintf('%1$s/run%2$04d_on.h5', cubedir, run);
    d.on.scan_var = h5read(fname_on, '/scan_var');
    d.on.imgs = h5read(fname_on, '/imgs');
    d.on.i0 = h5read(fname_on, '/i0');
    d.on.bin_counts = h5read(fname_on, '/bin_counts');    

    fname_off = sprintf('%1$s/run%2$04d_off.h5', cubedir, run);
    d.off.scan_var = h5read(fname_off, '/scan_var');
    d.off.imgs = h5read(fname_off, '/imgs');
    d.off.i0 = h5read(fname_off, '/i0');
    d.off.bin_counts = h5read(fname_off, '/bin_counts');    
end