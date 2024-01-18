function [equivvecs, min_distances] = findinwscell(r, A_s, prec)
% finds the equivalent of `r` in the Wigner-Seitz cell of columns of `A_s`.
% vectors with norms^2 within `prec` are considered equivalent.

%% all neighbors
nn =[-1    -1    -1
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

coordr=A_s\r;
allcoor = bsxfun(@plus, coordr - round(coordr), nn);
%
R = A_s*allcoor;
allnorms = (sum(R.^2,1));
%
[dmin, imin] = min(allnorms);
idxkeep=abs(dmin - allnorms) < prec;
equivvecs = allcoor(:,idxkeep);
min_distances = sqrt(allnorms(idxkeep));

end