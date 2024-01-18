function [vertices, polygonindices] = generate_face_vertices(A,b)
% FIXME: should remove duplicates out of vertices first, then crop the same
% elements out of polygonindices as well...

%to find the BZ vertices:
%for each set of three lattice vectors, look for the intersection point
%between the three normal planes.

eps_ = 1e-6;
[m,n]=size(A);
vertices=[];
polygonindices=[];
% vertices=zeros(n,1);
% polygonindices=zeros(n,1);
for i=1:(m-2)
    for j=(i+1):(m-1)
        for k=(j+1):m
            %find the point defined by A[i,j,k]*x == b[i,j,k]
            %if there is a solution, this is a vertex of the BZ
            try
                % each i, j or k represents a boundary face. The
                % intersection is the vertices
                x=A([i j k],:)\b([i j k]);
                if and(min((A*x-b))>-eps_,min((A*x-b))<Inf)
                    vertices=[vertices, x];
                    polygonindices=[polygonindices, [i, j, k]'];
                end
            end
        end
    end
end
vertices = vertices';
polygonindices=polygonindices';
% remove duplicate vertices
% [vertices, ia,ic] = unique(vertices,'rows');
% [vertices, ia] = remove_duplicates(vertices, 1e-6);
% polygonindices = polygonindices(ia,:);
end

