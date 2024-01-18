function psps = extract_lineout(imgs, cx, cy)
numbins = size(imgs,3);
psps = zeros(length(cx),numbins);
for ii = 1:size(imgs,3)
    tmp = imgs(:,:,ii);
    tmp = interp2(1:size(tmp,2),1:size(tmp,1),tmp,cx,cy);
    psps(:,ii) = tmp;
end
psps = psps';
end
