function allhkl = findclosesthkl2(bg,Q)
%% define nearest neighbors
nnhkls =[-1    -1    -1
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
nng = nnhkls*bg';

%% find nearest neighbor to `Q` by checking next neighbors defined in `nns`
qrd  = round(Q*inv(bg)')*bg';
allhkl = zeros(size(Q));
for k=1:size(Q,1)
    nn = bsxfun(@plus, qrd(k,:), nng);
    R2 = sum(nn.^2,2)';
    QdR= 2*Q(k,:)*nn';
    dis = bsxfun(@minus, R2, QdR);
    [D, ind] = min(dis,[],2);
    allhkl(k,:) = nn(ind,:);
end
allhkl = round(allhkl*inv(bg)');
1;
