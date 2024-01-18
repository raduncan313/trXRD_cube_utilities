[hkls Ks] = plotdata(Q,target_fun,allhkl,Rot,[0 2]);
%%
%figure
hold on
for ii=1:size(Ks,2)
	scatter3(Ks(1,ii),Ks(2,ii),Ks(3,ii),75,'filled','red');
end

%%
set(gca,'nextplot','replacechildren');
set(gca,'DataAspectRatio',[1 1 1]); % 'PlotBoxAspectRatio',[1 1 1],'ZLim',[0 0.8]);
set(gca,'CameraViewAngleMode','manual');
view(0,15);

%%
P = [0 0 1];
verts = [];
count = 0;
for ii = 1:size(Ks,2)
	count = count +1;
	verts(:,count) = Ks(:,ii)+Rot*P'/a0;
end

hold on
for ii=1:size(verts,2)
	scatter3(verts(1,ii),verts(2,ii),verts(3,ii),75,'filled','blue');
end
