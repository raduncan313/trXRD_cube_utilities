function [uA, id] = remove_duplicates(A, eps_)
% function [uA, id] = remove_duplicates(A, eps_)
% assumes columns of A are n-vectors and removes duplicate columns
% eps_ = 1e-6;  % precision
[m,n]=size(A);
% At = A';

del = zeros(m,m);

for i=1:m
    for j = i+1:m
        dif = A(i,:) - A(j,:);
        if norm(dif) < eps_
            del(i,j) = 1;
        end
    end
end

% del = bsxfun(@minus, A, At); %...
id = sum(del,2)==0;
uA = A(id,:);

end
