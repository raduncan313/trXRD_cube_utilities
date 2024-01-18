function [p, I, d] = nearest_neighbors(P, S, varargin)
% function [id, d] = NEAREST_NEIGHBORS(P, S)
%
% finds the nearest neighbors of P in the set S. 
%
% Input:
%   P  D x N vector of points to query
%   S  D x NS vector of points in the set to search
% Output:
%   p  list of nn's found 
%   I  indices of p in the columns of S

% is_number=true;
% eps_ = 1e-6;

% [num_nn] = check_inputs(varargin);

[D, N] = size(P);
[D2,NS] = size(S);
if D ~= D2
    error('Nearest_neighbors: Input dimensions must match.')
end
p=[];
I=[];
d=[];
% if is_number
    for i=1:N
        dis = zeros(NS,1);
        for j=1:NS
            dis(j) = sum((P(:,i) - S(:,j)).^2,1);
        end
        [s_dis, id] = sort(dis, 'ascend');
        s_dis = s_dis(1:num_nn);
%         id = id(1:num_nn);
        idx=id(1:num_nn);
        p = [p, S(:,id)];
        I = [I, idx'];
        d = [d, s_dis'];
    end
% else
%     for i=1:N
%         dis = zeros(NS,1);
%         for j=1:NS
%             dis(j) = sum((P(:,i) - S(:,j)).^2,1);
%         end
%         min_d = min(dis);
%         id = abs(dis-min_d) < eps_;
%         idx = find(id);
%         p = [p, S(:,idx)];
%         I = [I, idx'];
%         d = [d, dis(idx)'];
%     end
end

function num_nn = check_inputs(vararg)
    % defaults:
    if nargin == 1
%         is_number=true;
%         epsilon = 1e-6; % not used when `is_number == true`
        num_nn=vararg{1};
    end
    
end
        