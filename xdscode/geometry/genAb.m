function [A, b] = genAb(bg, nr1, nr2, nr3)
% `genAb(bg, nr1 = 1, nr2 = 1, nr3 = 1)`

% Generates a list of neighbors `A` from a basis given by the columns of
% `bg`. Also returns the norm^2/2 of each vector in `b`. Useful for
% determining the Wigner-Seitz cell.

% # nr1,nr2,nr3=1,1,1 #
A=zeros((2*nr1+1)*(2*nr2+1)*(2*nr3+1),3);
b=zeros((2*nr1+1)*(2*nr2+1)*(2*nr3+1),1);
ct=1;

for n1=-nr1:nr1
    for n2=-nr2:nr2
        for n3=-nr3:nr3
            b(ct) = -norm(bg*[n1; n2; n3]).^2/2;
            A(ct, :) = bg*[n1; n2; n3];
            ct = ct + 1;
        end
    end
end
%%
% [d, id] = sort(b);
% A = A(id,:);
