function [Q, QQ, allK_q, allhkl] = generate_reduced_wavevectors(geometry)
%%
lambda0 = geometry.lambda0;
ny=geometry.imageNy;
nz=geometry.imageNz;

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
if ~isfield(geometry,'rot_matrix')
    geometry.rot_matrix = @huber_matrix;
end
rot_matrix = geometry.rot_matrix;
Rot = rot_matrix(phi0, theta, chi)*SamRot;

bg = geometry.primvects';

% Q = det_kspace_proj(geometry);
Q = det_kspace_sphere(geometry);
Q=reshape(Q,[],3);

%% generate a longer list of HKLs candidate, will restrict later
% unihkls = unique(round(inv(bg)*inv(Rot)*Q')','rows');  % unique HKL indices 
% gunihkls=unihkls*bg';                    % Cartesian coords. of unihkls
% morehkls = [];
% nn=3;
% for x=-nn:nn
%     for y = -nn:nn
%         for z=-nn:nn
%             for k=1:size(unihkls,1)
%                 morehkls = [morehkls; [x,y,z] + unihkls(k,:)];
%             end
%         end
%     end
% end
% morehkls = unique(morehkls, 'rows');
% gmorehkls = morehkls*bg';

%% plot HKL candidates
% figure(16); clf
% plot_tds_on_surface(reshape(Q*Rot, nz,ny,3),zeros(nz,ny)); 
% hold all
% scatter3(gunihkls(:,1),gunihkls(:,2),gunihkls(:,3),50,'b','filled');
% scatter3(gmorehkls(:,1),gmorehkls(:,2),gmorehkls(:,3),'r')

%% remove absolutely unreachable HKLs
% loop over  a bunch of rec.lat. points close the the Ewald sphere 
% 
% k0 = 1/lambda0*[1, 0, 0]*inv(Rot)'; % in same coordinates as `at` and `bg`
% bmin = max(sum(bg.^2,1));
% k00 = 1/lambda0;
% 
% gnnhkls = [];
% for k = 1:size(gmorehkls,1)
%     thisg = gmorehkls(k,:);  % replace with `gmorehkls` if more neighbors needed
%     kp = thisg+k0;
% %     if (	(abs(norm(kp) - k00) < 5*bmin)  )
%         gnnhkls = [gnnhkls; thisg];
% %     end
% end
% nnhkls = (bg\gnnhkls')';

%% plot extended list of candidates
% figure(17);
% plot_tds_on_surface(reshape(Q*Rot, nz,ny,3),zeros(nz,ny)); 
% hold all;
% scatter3(gnnhkls(:,1),gnnhkls(:,2),gnnhkls(:,3),'k')

%% A few ways to find nearest neighbors.  
% for each Q point, searches which of verts(:) is closest then saves the
% indices. Then ind(i) is the index of verts closest to Q(:,i)
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by eahc row of Q (m x d array).

% [ind, D] = knnsearchFEX(Q,gnnhkls);	

%% search for nearest rec.lat vector in the set `gnnhkls` for each Q(i,:)
% nq = size(Q,1);
% ind = zeros(nq,1);
% Q2 = sum(Q.^2,2);
% R2 = sum(gnnhkls.^2,2);
% QdR= Q*gnnhkls';
% for k=1:nq;
%     currdist = (Q(k,1)-gnnhkls(:,1)).^2 + (Q(k,2)-gnnhkls(:,2)).^2 + (Q(k,3)-gnnhkls(:,3)).^2;
% %     currdist = Q2(k) + R2 - QdR(k,:)';
%     [D, indmin] = min(currdist);
%     ind(k)= indmin;
% end

%% faster way to search for NN
% nq = size(Q,1);
% % K = size(gnnhkls,1);
% % Q2 = sum(Q.^2,2);
% R2 = sum(gnnhkls.^2,2)';
% % [~,ind] = sort(R2);
% QdR= 2*Q*gnnhkls';
% % dis = Q2(:,ones(1, K)) + R2(ones(nq, 1),:) - QdR; % <-- the minimum of this
% % dis = R2(ones(nq, 1),:) - QdR;                    % <-- is the same as the min of this
% dis = bsxfun(@minus, R2, QdR);
% [D, ind] = min(dis,[],2);

%% simple is usually faster:
% Qrot = Q*inv(Rot)'*inv(bg)';
allhkl = findclosesthkl2(bg, Q*inv(Rot)');
allK_q = allhkl*(bg');
% profile viewer

%%
%save the closest rec.lat. vector to each pixel
% allK_q = gnnhkls(ind,:);
%find the reduced wavevector
QQ = (Q*inv(Rot)' - allK_q);            % Reduced wavevector [A^-1]
% allhkl = nnhkls(ind,:);				% also the hkl indices
1;
