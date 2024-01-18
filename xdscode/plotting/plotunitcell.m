function bz=plotunitcell(bg, varargin)
    dim = size(bg,1);

    if nargin<2
        t0 = zeros([dim,1]);
    else
        t0 =varargin{1};
    end
	ibg = inv(bg)';
    for i=1:dim 
        ibg(:,i)= ibg(:,i)/norm(ibg(:,i));
    end
% 	b= vcat(-sum(ibg.*(bg.+t0),1)[:], sum(ibg.*t0,1)[:])
    v1=-sum(ibg.*bsxfun(@plus,bg,t0),1)';
    v2=sum(bsxfun(@times,ibg, t0),1)';
    b = cat(1,v1,v2);
% 	b= cat(-sum(ibg.*(bg.+t0),1), sum(ibg.*t0,1));
	A= cat(1,-ibg', ibg');
	patches=plotregion(A,b,[],[],[0.9 0.9 0.9],0.2);
%     axis equal tight vis3d
%     hold on

bz=hggroup(gca);
set(findobj(gca,'type','line'),'parent',bz); % put all line plots into group
for ii=1:length(patches)
    set(patches(ii),'Parent',bz)
end

axis tight equal vis3d

rotate3d;

end
