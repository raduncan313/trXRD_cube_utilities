function d = norm_i0(d)
    sz = size(d.on.imgs);
    d.on.imgs = d.on.imgs./permute(repmat(d.on.i0, [1, sz(1), sz(2)]), [2,3,1]);
    d.off.imgs = d.off.imgs./permute(repmat(d.off.i0, [1, sz(1), sz(2)]), [2,3,1]);    
end