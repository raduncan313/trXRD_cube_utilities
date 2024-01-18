function [ww, vv] = phonon_freqs_vecs(Dmats)
nq = size(Dmats, 3);
nd = size(Dmats,1);
%%
ww=zeros(nd, nq);
vv=zeros(nd,nd,nq);
for ii=1:nq
    DD = Dmats(:,:,ii);
    [vv(:,:,ii), w2] = eig((DD+DD')/2);
    w2 = diag(w2);
    negid = (w2 < 0);
    w2(negid) = -w2(negid);
    w = sqrt(w2);
    w(negid) = -w(negid);
    ww(:,ii) = w;
end

end
