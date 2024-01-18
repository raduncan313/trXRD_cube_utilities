function plotbzmesh(geometry,data)
geometry.imageNz = floor(size(data,1)/4);
geometry.imageNy = floor(size(data,2)/4);
gnz = geometry.imageNz;
gny = geometry.imageNy;
nz = size(data,1);
ny = size(data,2);

[Q, QQ, ~, allhkl] = generate_reduced_wavevectors(geometry);
hold on; 
allhkl = reshape(allhkl,[gnz, gny, 3]);
% Q=reshape(Q,[gnz, gny, 3]);
% QQ=reshape(QQ,[gnz, gny, 3]);
rr = imresize(abs(sqrt(2)*allhkl(:,:,1) - pi*allhkl(:,:,2) + exp(1)*allhkl(:,:,3)), [nz, ny]);
% rr=abs(reshape(sqrt(2)*allhkl(:,1) - pi*allhkl(:,2) + exp(1)*allhkl(:,3), [Nz, Ny]));
c1=contour(rr,100,'color','b','linewidth',2);
hold off
% imagesc(rr);
% 
% allhkl=reshape(allhkl,[Nz, Ny, 3]);
% lambda0=geometry.lambda0;
% bg = geometry.primvects';
% 
% theta   = geometry.theta;
% chi     = geometry.chi;
% phi0    = geometry.phi;
% SamRot 	= geometry.SamRot;
% Rot = huber_matrix(phi0, theta, chi)*SamRot;

%% set the datacursor to update with useful information
% dcm_obj = datacursormode(gcf);
% set(dcm_obj, 'DisplayStyle', 'window')
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,{QQ,Q,allhkl},lambda0,bg,Rot})

1;

function txt = myupdatefcn(~,event_obj,arrbundle,lambda0,bg,Rot)  %
QQ = arrbundle{1};
Q = arrbundle{2};
allhkl = arrbundle{3};
% Customizes text of data tips
pos = get(event_obj,'Position');
% I = get(event_obj, 'DataIndex');
q = squeeze(Q(pos(2), pos(1), :));
qq= squeeze(QQ(pos(2), pos(1), :));
hkl=squeeze(allhkl(pos(2), pos(1), :));
MM = inv(bg)*inv(Rot);
k0 = [1,0,0]'/lambda0;
kp = k0+q(:);
mu = atan2d(kp(2), kp(1));
del = 90-acosd(kp(3)*lambda0);
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       sprintf('hkl: %i, %i, %i',hkl),...
       sprintf('mu: %g    del: %g', mu, del),...  %['mu: ', num2str(mu), '\tdel: ', num2str(del)]
       sprintf('ThetaB: %g    inv(d): %g', asind(lambda0*norm(q)/2), 1/norm(q)),...
       sprintf('q : %1.3g  %1.3g  %1.3g', MM*qq),...
       sprintf('Q : %1.3g  %1.3g  %1.3g', inv(bg)*q)};
       