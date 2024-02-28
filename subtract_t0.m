function d = subtract_t0(d, t0)
    d.on.scan_var = d.on.scan_var - t0;
    d.off.scan_var = d.off.scan_var - t0;
end