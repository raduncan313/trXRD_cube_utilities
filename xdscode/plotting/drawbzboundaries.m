function [dda_idy, dda_idx]=drawbzboundaries(hkls, m, n)
    aa=reshape([1 1/sqrt(2) 1/sqrt(3)]*hkls', m, n);
    da1=diff(aa,1);
    da2=diff(aa,[],2);
    dda=zeros(size(aa));
    dda(1:end-1,:) = abs(da1);
    dda(:,1:end-1) = dda(:,1:end-1) + abs(da2);
    [dda_idx, dda_idy] = find(dda);
    plot(dda_idy, dda_idx, 'g.')
end
