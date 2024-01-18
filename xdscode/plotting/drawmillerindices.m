function uhkls=drawmillerindices(hkls, nqz, nqy)
    uhkls=unique(hkls,'rows');
    for k=1:size(uhkls,1)
        aa=reshape( all(uhkls(k,:) == hkls,2), nqz, nqy);
        [idr, idc]=find(aa);
        cmin = min(idc);
        cmax = max(idc);
        rmin = min(idr);
        rmax = max(idr);
        text(cmin+(cmax-cmin)/3, rmin+(rmax-rmin)/2, num2str(uhkls(k,:)),...
            'color','g','fontsize',12);
    end
end
