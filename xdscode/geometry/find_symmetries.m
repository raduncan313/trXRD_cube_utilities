function rotations = find_symmetries(at, prec)
% function rotations = find_symmetries(at)
% Finds all allowed proper and improper rotations allowed by the symmetry
% of the unit cell, specified as the columns of `at`. Comparing the rotated
% structure is done within the precision specified by `prec`
% 

if nargin < 2
    prec = 1e-6;
end
    
l=1;
    relative_axes = zeros(3, 27);
    for i=1:3
        for j=1:3 
            for k=1:3
               relative_axes(:,l) = [i,j,k] - 2;
               l = l + 1;
            end
        end
    end

    orig_metric = at'*at;
    
	ct = 0;
	numsyms=0;
	rotations = [];
	for k=1:27
        for j=1:27
            for i=1:27
        		axes = [relative_axes(:, i) relative_axes(:, j) relative_axes(:, k)];
                
                if abs(det(axes)) == 1
        			lattice = at*axes;
        			new_metric = lattice'*lattice;
%         			if (is_same_metric(new_metric, orig_metric, 1e-6))
                    if norm(new_metric-orig_metric) < prec
                        rotations = [rotations, axes];
        				numsyms=numsyms+1;
                    end
        			ct = ct + 1;
                end
            end
        end
    end
    rotations = reshape(rotations, 3,3,[]);
    
end


function tf = is_same_metric(Gorig, Gnew, prec)
% 	# since G_{i,j} = a_i^T * a_j :
% 	# cos(anglebetween(a_i, a_j)) = G_ij /sqrt(G_{i,i})/sqrt(G_{j,j}) 
    if ~all(abs(diag(Gorig) - diag(Gnew)) < prec)
    	tf = false;
        return;
    end

	for i=1:2
		for j=(i+1):3
            % compares the off-diagonal 
            tf = (abs(Gorig(i,j) - Gnew(i,j)) < prec);
%             tf = (abs(get_angle(Gorig,i,j) - get_angle(Gnew,i,j)) < prec);
            if ~tf
                return
            end
			% Spglib uses the delta_theta between the two angles
			% then if sin(dtheta)^2 != 0
			% compares sin_dtheta2 * length_ave2 > symprec^2
			% (what is symprec???? has units??)
		end
    end
    norm(Gorig-Gnew);
    tf = true;
end

% function ang = get_angle(metric, i, j)
%     n1  = sqrt(metric(i,i));
%     n2  = sqrt(metric(j,j));
%     ang = acos(metric(i,j)/n1/n2);
% end
% 
