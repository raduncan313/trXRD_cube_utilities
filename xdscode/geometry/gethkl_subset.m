function varargout = gethkl_subset(allhkl,hkl)
% returns the indices of the 2D image and the reduced wavevectors that are
% closest to a certain [h k l]

Kmask = squeeze(allhkl(1,:,:)==hkl(1) & allhkl(2,:,:)==hkl(2) & allhkl(3,:,:)==hkl(3));
%QQr = QQ(:,Kmask);
%scatter3(QQ(1,Kmask),QQ(2,Kmask),QQ(3,Kmask),0.2,'filled')
if nargout == 0
%	scatter3(a0*QQ(1,Kmask),a0*QQ(2,Kmask),a0*QQ(3,Kmask),0.2,'filled')
	imagesc(reshape(Kmask,[512, 512]) ); axis image
else
	varargout{1} = Kmask;
end
