function [Q, allhkl] = plotdata(geometry,data,varargin)
% some help here

geometry.imageNz = size(data,1);
geometry.imageNy= size(data,2);

Nz = geometry.imageNz;
Ny = geometry.imageNy;

% plots the ratio data on the detector scale (pixels):
if isempty(varargin)
    imagesc(data); axis image;
else
    imagesc(data,varargin{1}); axis image;
end

[Q, QQ, allK_q, allhkl] = generate_reduced_wavevectors(geometry);
Q=reshape(Q,[Nz, Ny, 3]);
QQ=reshape(QQ,[Nz, Ny, 3]);
allK_q=reshape(allK_q,[Nz,Ny,3]);

hold on

% rr=abs(reshape(sqrt(2)*allhkl(:,1) - pi*allhkl(:,2) + exp(1)*allhkl(:,3), [Nz, Ny]));
% contour(rr,30,'color','w');

drawbzboundaries(allhkl, Nz, Ny);
drawmillerindices(allhkl, Nz, Ny);

if size(varargin,2) > 1
    axis(varargin{2})
else
    axis image
end

allhkl=reshape(allhkl,[Nz, Ny, 3]);
k0 = [1,0,0]'/geometry.lambda0;
lambda0=geometry.lambda0;
bg = geometry.primvects';

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
if ~isfield(geometry,'rot_matrix')
    geometry.rot_matrix = @huber_matrix;
end
rot_matrix = geometry.rot_matrix;
Rot = rot_matrix(phi0, theta, chi)*SamRot;
MM = inv(bg)*inv(Rot);

%% set the datacursor to update with useful information
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'DisplayStyle', 'window')
set(dcm_obj,'UpdateFcn',{@myupdatefcn,{QQ,Q,allhkl},geometry,bg,Rot})

%% use ginput to set BZ labels on graph
%{
while (xi < Ny && yj < Nz)
	count = count +1;
%	figure(datafig);
	[x,y]=ginput(1);
	xi = round(x);	
	yj = round(y);
	try
		hkl = squeeze(allhkl(yj,xi,:));
        qq=squeeze(QQ(yj,xi,:));
        q=squeeze(Q(yj,xi,:));
        kp = k0+q;
        mu = atan2d(kp(2), kp(1));
        del = 90-acosd(kp(3)*geometry.lambda0);
 		text(x,y,sprintf('(%i %i %i)',hkl),'FontWeight','bold','Color','w','fontsize',14);
        
        disp(['(i, j) = (', num2str([xi, yj]),')' ])
        disp([sprintf('(%i %i %i)\t',hkl), ...
            'ThetaB = ', num2str(asind(lambda0*norm(q)/2)),...
            '    mu = ', num2str( mu), ...
            '    del = ', num2str(del) ])
        disp(['q = [',num2str((bg\qq)'), ']'])
        disp(['Q = [',num2str((MM*q)'), ']'])
        
    catch
       disp('Done.'); 
    end
end
%}
1;

function txt = myupdatefcn(~,event_obj,arrbundle,geometry,bg,Rot)  %
det_size_x = geometry.detector.det_pixels_horz;
det_size_y = geometry.detector.det_pixels_vert;
img_size_x = geometry.imageNy;
img_size_y = geometry.imageNz;

QQ = arrbundle{1};
Q = arrbundle{2};
allhkl = arrbundle{3};
lambda0 = geometry.lambda0;
theta = geometry.theta;
% Customizes text of data tips
pos = get(event_obj,'Position');
cdat=event_obj.Target.CData(pos(2),pos(1));
% I = get(event_obj, 'DataIndex');
q = squeeze(Q(pos(2), pos(1), :));
qq= squeeze(QQ(pos(2), pos(1), :));
hkl=squeeze(allhkl(pos(2), pos(1), :));
try
%     find_angles_outgoingk(geometry, hkl, theta);
end
% MM = inv(bg); %*inv(Rot)
MbtoB =([
    -1     1     1
     1    -1     1
     1     1    -1]);

k0 = [1,0,0]'/lambda0;
kp = k0+q(:);
mu = atan2d(kp(2), kp(1));
del = 90-acosd(kp(3)*lambda0);
txt = {['x (det X): ', num2str(pos(1)), '    (', num2str(pos(1)*det_size_x/img_size_x), ')'],...
       ['y (det Y): ', num2str(pos(2)), '    (', num2str(pos(2)*det_size_y/img_size_y), ')'],...
       sprintf('hkl: %i, %i, %i',hkl),...
       sprintf('mu: %g    del: %g', mu, del),...  %['mu: ', num2str(mu), '\tdel: ', num2str(del)]
       sprintf('2ThetaB: %g    d: %g', 2*asind(lambda0*norm(q)/2), 1/norm(q)),...
       sprintf('q : %1.3g  %1.3g  %1.3g', bg\qq),...
       sprintf('Q : %1.3g  %1.3g  %1.3g', (Rot*bg)\q),...
       sprintf('CData : %2.2E', cdat) };
%        sprintf('qc: %1.3g  %1.3g  %1.3g', MbtoB*(bg\qq)),...
       