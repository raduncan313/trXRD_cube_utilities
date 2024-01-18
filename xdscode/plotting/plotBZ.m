function bz=plotBZ(bg)

v1 = bg(:,1);
v2 = bg(:,2);
v3 = bg(:,3);

%
if size(v1,2)==1
	v1=v1';
end
if size(v2,2)==1
	v2=v2';
end
if size(v3,2)==1
	v3=v3';
end
%}

count = 0;
clear verts
%
NN=1;
for x= -NN:NN
	for y= -NN:NN
		for z= -NN:NN
			count = count +1;
			verts(count,:) = x*v1+y*v2+z*v3;
		end
	end
end

x0s = [];
for i=1:size(verts,1)
	x0s   = [x0s, norm(verts(i,:))];
end

persistent bznum;
if isempty(bznum)
    bznum=1;
end

% bzcolor = [0.9 0.2 0.1];
bzcolor=[0.85, 0.85, 0.85];

patches=plotregion(verts, -0.5*x0s.^2,[],[],bzcolor,0.1,[[0 0 0];v1;[0 0 0];v2;[0 0 0];v3]','r-');
axis equal tight; box on
set(gca,'CameraViewAngleMode','manual')
bz=hggroup(gca);
bz.DisplayName=['BZ',num2str(bznum)];
bznum=bznum+1;
set(findobj(gca,'type','line'),'parent',bz); % put all line plots into group
for ii=1:length(patches)
    set(patches(ii),'Parent',bz)
end
