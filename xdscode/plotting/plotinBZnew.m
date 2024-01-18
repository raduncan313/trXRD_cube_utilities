function allgroup=plotinBZnew(geometry,im0,hkl,varargin)
% plots the diffuse intensity in Intens inside the Brillouin zone for hkls
% for a given rotation matrix Rot (rotates b1,b2,b3). The function needs
% 'allhkl' to figure out which pixels of the image Intens belong to each 'hkls'
% zones.
nz = size(im0,1);
ny = size(im0,2);
geometry.imageNy=ny;
geometry.imageNz=nz;
hkl = hkl(:);
bg = geometry.primvects';
SamRot 	= geometry.SamRot;
rot_matrix = geometry.rot_matrix;

theta   = geometry.theta;
chi     = geometry.chi;
phi    = geometry.phi;
Rot = inv(rot_matrix(phi, theta, chi)*SamRot);

[Q, QQ, ~, allhkl] = generate_reduced_wavevectors(geometry);
% QQ = reshape(QQ, nz,ny,3);
Qr = reshape(Q, nz,ny,3);

% Q = det_kspace_proj(geometry);
% Q = reshape(Q,[],3);
% Qrel = Q/(bg/Rot);
% allhkl = findclosesthkl2(Qrel);
% Q = Q';
% QQ=QQ';
% allhkl=allhkl';
% QQ = reshape(bsxfun(@minus, (Q*Rot'), hkl'*bg), nz,ny,3);
% QQr = bsxfun(@minus, QQr, bg*hkl');

mask = squeeze(allhkl(:,1)==hkl(1) & allhkl(:,2)==hkl(2) & allhkl(:,3)==hkl(3));
mask = reshape(mask, nz,ny);
[rw,cl] = find(mask);
i1=min(rw);
i2=max(rw);
j1=min(cl);
j2=max(cl);

%%
QQr = Qr(i1:i2, j1:j2,:);
% QQr = QQ;
% im1=im0;
mask = mask(i1:i2, j1:j2,:);
im1 = im0(i1:i2, j1:j2,:);
im1(~mask) = nan;  % to ignore points outsize of the BZ domain

X = reshape(QQr(:,:,1),[],1);
Y = reshape(QQr(:,:,2),[],1);
Z = reshape(QQr(:,:,3),[],1);

%%
% X = reshape(Q(:,1),[],1);
% Y = reshape(Q(:,2),[],1);
% Z = reshape(Q(:,3),[],1);


%%
% QQ=reshape(QQ,3,[]);

%figure;
% clf
% looks for the BZ object to avoid replotting it
bz1=findobj(gca,'DisplayName',['BZ1',num2str(get(gcf,'Number'))]);
if isempty(bz1)
    bz1=plotBZ(bg);
    hold on
end

%% 

% generate evenly spaced points in the y-z plane
[yg, zg] = meshgrid(linspace(min(Y),max(Y),size(im1,1)),linspace(min(Z),max(Z),size(im1,2)));
% [xg, yg] = meshgrid(linspace(min(X),max(X),size(im1,1)),linspace(min(Y),max(Y),size(im1,2)));

% interpolate the QQ 'scattered' points onto a the evenly spaced grid
xg = griddata(QQr(:,:,2),QQr(:,:,3),QQr(:,:,1),yg,zg);
% zg = griddata(QQr(:,:,1),QQr(:,:,2),QQr(:,:,3),xg,yg);

% interpolate the intensity as if it were a function im0 = f(qy, qz)
cgrid = griddata(QQr(:,:,2),QQr(:,:,3),im1,yg,zg); %,'cubic'
% cgrid = griddata(QQr(:,:,2),QQr(:,:,3),im1,xg,yg); %,'cubic'


%% shift and rotate Q's to first BZ:
Ghkl = bg*hkl;
% allq = bsxfun(@minus, [xg(:), yg(:), zg(:)]*Rot', Ghkl');
allq = [xg(:), yg(:), zg(:)]*Rot';
xg = reshape(allq(:,1), size(QQr,2), size(QQr,1)) - Ghkl(1);
yg = reshape(allq(:,2), size(QQr,2), size(QQr,1)) - Ghkl(2);
zg = reshape(allq(:,3), size(QQr,2), size(QQr,1)) - Ghkl(3);
% % xg = reshape(allq(:,1) - Ghkl(1);
% % yg = yg - Ghkl(2);
% % zg = zg - Ghkl(3);

%% scratch
% nonnans = ones(size(zg));
% nonnans(isnan(xg)) = 0;
% nonnans(isnan(yg)) = 0;
% nonnans(isnan(zg)) = 0;
% nonnans(isnan(cgrid)) = 0;
% tri = delaunay(yg(nonnans>0),zg(nonnans>0));
% sanity check:
% trisurf(tri, xg(nonnans>0), yg(nonnans>0), zg(nonnans>0), cgrid(nonnans>0),'edgecolor','none');

%%

%%
% figure(5);clf
% thisSurf = surf(xgrid,yg,zg,(cgrid),'EdgeColor','none',...
thisSurf = surf(xg,yg,zg,(cgrid),'EdgeColor','none',...
								'FaceColor','interp',...
								'FaceAlpha',0.9,... 
								'FaceLighting','phong' );
%                                 'FaceLighting','gouraud');
%
axis tight vis3d ;
%%
% set(gca,'nextplot','replacechildren');
set(gca,'DataAspectRatio',[1 1 1]); % 'PlotBoxAspectRatio',[1 1 1],'ZLim',[0 0.8]);
set(gca,'CameraViewAngleMode','manual');
% view(0,15);
Gu = bg*hkl;
Gu = det(bg)^(1/3)*Gu/norm(Gu); % shrink to plot
quiv1=quiver3(0,0,0,Gu(1),Gu(2),Gu(3),'linewidth',1);
txt1=text(Gu(1),Gu(2),Gu(3), num2str(hkl'),'Color','k','fontsize',12);
% text(q(1,1),q(2,1),q(3,1), num2str(hkl'),'FontWeight','bold','Color','k','fontsize',14);

allgroup=hggroup(gca);
bz1.Parent=allgroup;
quiv1.Parent=allgroup;
txt1.Parent=allgroup;
thisSurf.Parent=allgroup;

%%
if nargin > 3
    caxis(varargin{1});
end
rotate3d on
%axis off
% colorbar([0.11 0.95 0.75 0.03])
