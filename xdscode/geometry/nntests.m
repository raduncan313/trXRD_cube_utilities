%% define nearest neighbor indices
nns =[-1    -1    -1
    -1    -1     0
    -1    -1     1
    -1     0    -1
    -1     0     0
    -1     0     1
    -1     1    -1
    -1     1     0
    -1     1     1
     0    -1    -1
     0    -1     0
     0    -1     1
     0     0    -1
     0     0     0
     0     0     1
     0     1    -1
     0     1     0
     0     1     1
     1    -1    -1
     1    -1     0
     1    -1     1
     1     0    -1
     1     0     0
     1     0     1
     1     1    -1
     1     1     0
     1     1     1];


%%
Q = det_kspace_proj(geometry); 
sizq=size(Q);
bg = geometry.primvects';
% Q = reshape(Q,[],3);

%%
Q = det_kspace_proj(geometry); 
Q = reshape(Q,[],3)*inv(Rot)';
qrd  = round(Q*inv(bg)');
allhkl = zeros(size(Q));
allKq = zeros(size(Q));
tic
% profile on
for k=1:size(Q,1);
%     [allKq(k,:), allhkl(k,:)]=findclosesthkl(bg, Q(k,:));
%     allhkl(k,:)=findclosesthkl(bg, Q(k,:));
    hkl = bsxfun(@plus, qrd(k,:), nns)*bg';
    R2 = sum(hkl.^2,2)';
    QdR= 2*Q(k,:)*hkl';
    dis = bsxfun(@minus, R2, QdR);
    [D, ind] = min(dis,[],2);
    allhkl(k,:) = hkl(ind,:);
end
% profile viewer
toc

%% test cache locality with transposed vectors
nnst =[-1    -1    -1
    -1    -1     0
    -1    -1     1
    -1     0    -1
    -1     0     0
    -1     0     1
    -1     1    -1
    -1     1     0
    -1     1     1
     0    -1    -1
     0    -1     0
     0    -1     1
     0     0    -1
     0     0     0
     0     0     1
     0     1    -1
     0     1     0
     0     1     1
     1    -1    -1
     1    -1     0
     1    -1     1
     1     0    -1
     1     0     0
     1     0     1
     1     1    -1
     1     1     0
     1     1     1]';


Q = det_kspace_proj(geometry); 
Q = bg\reshape(Q,[],3)';
qrd  = round(Q);
allhkl = zeros(size(Q));
allKq = zeros(size(Q));
tic
% profile on
for k=1:size(Q,2);
    hkl = bsxfun(@plus, qrd(:,k), nnst);
    R2 = sum(nnst.^2,1);
    QdR= 2*Q(:,k)'*hkl;
    dis = bsxfun(@minus, R2, QdR);
    [D, ind] = min(dis,[],2);
    allhkl(:,k) = hkl(:,ind);
end
% profile viewer
toc

%%

theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;
bg=geometry.primvects'

Q = det_kspace_proj(geometry); 
Qrel = reshape(Q,[],3)/(bg/Rot);
% profile on
allhkl = findclosesthkl2(bg,Qrel);
% profile viewer


%%
profile on
bzs=makebzmesh(geometry, [nz,ny]);
profile viewer