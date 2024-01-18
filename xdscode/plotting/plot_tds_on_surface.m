function plot_tds_on_surface(Q,im0)

ny=size(Q,2);
nz=size(Q,1);
X = reshape(Q(:,:,1),[],1);
Y = reshape(Q(:,:,2),[],1);
Z = reshape(Q(:,:,3),[],1);

% creates a grid where to evaluate the interpolated surface

[yg, zg] = meshgrid(linspace(min(Y),max(Y),ny),linspace(min(Z),max(Z),nz));

% interpolates the surface defined by points [X Y Z] to the domain [xg yg],
% stores it into xgrid = f(yg,zg)
%xgrid = griddata(Y,Z,X,yg,zg);
xgrid = griddata(squeeze(Q(:,:,2)),squeeze(Q(:,:,3)),squeeze(Q(:,:,1)),yg,zg);

% also need to interpolate the image to the evenly spaced grid!!!!
cgrid = griddata(squeeze(Q(:,:,2)),squeeze(Q(:,:,3)),im0,yg,zg);


%% plots the interpolated image as caxis on the 3D surface
%surf(yg,zg,xgrid,flipud(Iss{1}),'EdgeColor','none'); caxis([0 1e3]);
%set(gca,'CameraViewAngleMode','auto');

%figure;
%clear Mov;
thisSurf = surf(xgrid,yg,zg,(cgrid),'EdgeColor','none',...
								'FaceColor','interp',...
								'FaceAlpha',0.75,... 
								'FaceLighting','phong' );

axis tight vis3d off;

% set(gca,'nextplot','replacechildren');
set(gca,'DataAspectRatio',[1 1 1]); % 'PlotBoxAspectRatio',[1 1 1],'ZLim',[0 0.8]);
set(gca,'CameraViewAngleMode','manual');
view(0,15);
