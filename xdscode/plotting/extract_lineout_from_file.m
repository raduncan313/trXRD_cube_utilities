function [delays, psps, pspsnorm, imgs, normimgs] = extract_lineout_from_file(filename,cx,cy,dataopts)
rundata = load(filename,'dd');
[imgs, delays] = normalize_i0(rundata.dd);
imgs = smoothdata(imgs, 1, 'gaussian',dataopts.img_smooth_pixel);
imgs = smoothdata(imgs, 3, 'gaussian',dataopts.img_smooth_time);
normimgs = imgs./mean(imgs(:,:,dataopts.meanrange),3);
psps = extract_lineout(imgs, cx, cy)';
pspsnorm = extract_lineout(normimgs, cx, cy)';
