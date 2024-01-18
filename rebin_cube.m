function d = rebin_cube(d0, binint)

    scan_var_on_0 = d0.on.scan_var;
    imgs_on_0 = d0.on.imgs;
    i0_on_0 = d0.on.i0;
    bin_counts_on_0 = d0.on.bin_counts;
    
    scan_var_off_0 = d0.off.scan_var;
    imgs_off_0 = d0.off.imgs;
    i0_off_0 = d0.off.i0;
    bin_counts_off_0 = d0.off.bin_counts;

    len1 = length(imgs_on_0(:,1,1));
    len2 = length(imgs_on_0(1,:,1));
    len3_0 = length(imgs_on_0(1,1,:));
    len3 = floor(len3_0/binint);

    scan_var_on = zeros(len3,1);    
    imgs_on = zeros(len1,len2,len3);
    i0_on = zeros(len3,1);
    bin_counts_on = zeros(len3,1);
    
    scan_var_off = zeros(len3,1);
    imgs_off = zeros(len1,len2,len3);
    i0_off = zeros(len3,1);
    bin_counts_off = zeros(len3,1);

    hh = 1;
    for ii = 1:len3
        imb_on = zeros(len1,len2);
        i0b_on = 0;
        bcb_on = 0;
        
        imb_off = zeros(len1,len2);
        i0b_off = 0;
        bcb_off = 0;
        
        for jj = 1:binint
            imb_on = imb_on + imgs_on_0(:,:,hh);
            i0b_on = i0b_on + i0_on_0(hh);
            bcb_on = bcb_on + bin_counts_on_0(hh);
            
            imb_off = imb_off + imgs_off_0(:,:,hh);
            i0b_off = i0b_off + i0_off_0(hh);
            bcb_off = bcb_off + bin_counts_off_0(hh);
            
            hh = hh + 1;
        end
        
        imgs_on(:,:,ii) = imb_on;
        i0_on(ii) = i0b_on;
        bin_counts_on(ii) = bcb_on;
        
        imgs_off(:,:,ii) = imb_off;
        i0_off(ii) = i0b_off;
        bin_counts_off(ii) = bcb_off;
        
        scan_var_on(ii) = mean([scan_var_on_0(hh-1), scan_var_on_0(hh-binint)]);
        scan_var_off(ii) = mean([scan_var_off_0(hh-1), scan_var_off_0(hh-binint)]);
    end
    
    d.on.scan_var = scan_var_on;
    d.on.imgs = imgs_on;
    d.on.i0 = i0_on;
    d.on.bin_counts = bin_counts_on;
    
    d.off.scan_var = scan_var_off;
    d.off.imgs = imgs_off;
    d.off.i0 = i0_off;
    d.off.bin_counts = bin_counts_off;
    
end