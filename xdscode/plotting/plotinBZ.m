function plotinBZ(geometry,Intens,hkls,varargin)
% plots the diffuse intensity in Intens inside the Brillouin zone for hkls
% for a given rotation matrix Rot (rotates b1,b2,b3). The function needs
% 'allhkl' to figure out which pixels of the image Intens belong to each 'hkls'
% zones.

geometry.imageNy=size(Intens,1);
geometry.imageNz=size(Intens,2);

%{
if size(b1,1)==1
	b1=b1';
end
if size(b2,1)==1
	b2=b2';
end
if size(b3,1)==1
	b3=b3';
end

b1 = Rot*b1;
b2 = Rot*b2;
b3 = Rot*b3;
%}
[~, QQ, ~, allhkl] = generate_reduced_wavevectors(geometry);
%Intens(isnan(Intens)) = 0;
% QQ=reshape(QQ,3,[]);

%figure;
% plotBZ(a0*Rot*b1,a0*Rot*b2,a0*Rot*b3);

bg = geometry.primvects';
SamRot 	= geometry.SamRot;

theta   = geometry.theta;
chi     = geometry.chi;
phi    = geometry.phi;
Rot = inv(huber_matrix(phi, theta, chi)*SamRot);
plotBZ(bg);

% hold on

for l=1:size(hkls,2)
	Hhkl = hkls(:,l);
%     mask = gethkl_subset(allhkl,Hhkl);
    mask = squeeze(allhkl(:,1)==hkls(1,l) & allhkl(:,2)==hkls(2,l) & allhkl(:,3)==hkls(3,l));

	q = QQ(mask,:);
	scatter3(q(:,1),q(:,2),q(:,3),50,Intens(mask'),'filled');
    Gu = bg*Hhkl;
    Gu = det(bg)^(1/3)*Gu/norm(Gu); % shrink to plot
    quiver3(0,0,0,Gu(1),Gu(2),Gu(3),'linewidth',1)
    text(Gu(1),Gu(2),Gu(3), num2str(Hhkl'),'FontWeight','bold','Color','k','fontsize',14)
    text(q(1,1),q(1,2),q(1,3), num2str(Hhkl'),'FontWeight','bold','Color','k','fontsize',14);
end
caxis(varargin{1});
%axis off
% colorbar([0.11 0.95 0.75 0.03])
