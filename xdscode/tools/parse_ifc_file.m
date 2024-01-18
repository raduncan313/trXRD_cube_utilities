function [at, tau, masses, nr1, nr2, nr3, phid] = parse_ifc_file(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parse full file using fopen and fread (you should use this version)
bohr_angstr = 0.52917720859;
amu_to_ry=911.4442;
fp=fopen(filename,'r'); % output of q2r.x

A=fscanf(fp,'%d %d %d %f %f %f %f %f %f',9);
ntyp = A(1);
nat =  A(2);
ibrav = A(3); % bravais lattice index

% the FCC unit cell in QE is defined as this:
%       v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)
a0 = A(4)*bohr_angstr;
if ibrav == 2
    at = 0.5*a0*[-1 0   1
                 0  1   1
                -1  1   0]';

else %FIXME: need to update if we ever use other systems
    at = eye(3)*a0;
end

% read types and names of atoms
for ii=1:ntyp
    idx=fscanf(fp,'%d',1);
    atm=fscanf(fp,'%s',2);
    amass(ii)=fscanf(fp,'%f',1);
    %TODO: check if idx = ii
end

% read atomic positions in the unit cell
for ii=1:nat
    A=fscanf(fp,'%d %d %f %f %f',5);
    idx = A(1);
    ityp(ii) = A(2);
    tau(:,ii) = A(3:end);
end
masses=amass(ityp)/amu_to_ry;

% check and read if has tensors
has_zstar=fscanf(fp,'%s',1);
if has_zstar=='T'
    epsil=fscanf(fp,'%f %f %f %f %f %f %f %f %f',[3,3]);
    for na=1:nat
        idx=fscanf(fp,'%d',1);
        zeu(:,:,na)=fscanf(fp,'%f %f %f %f %f %f %f %f %f',[3,3]);
    end
end

% number of neighbors
A=fscanf(fp,'%d %d %d',3);
nr1=A(1); nr2=A(2); nr3=A(3);

%************************** read IFCs from file 

si_ifc = fscanf(fp, '%f ',[4, (nr1*nr2*nr3+1)*3*3*nat*nat])';
%
phid = zeros(3,3,nat,nat,nr1,nr2,nr3);
mat_ind=zeros(3,3,nat,nat,4);
count = 1;
for i=1:3
    for j=1:3
        for nb=1:nat
            asr_sum=0;
            for na=1:nat
%                 count
                mat_ind(i,j,na,nb,:) = si_ifc(count,:);
                count = count+1;
                for n3=1:nr3
                    for n2=1:nr2
                        for n1=1:nr1
%                             atpos(n1,n2,n3) = si_ifc(count,1:3);
                            phid(i,j,na,nb,n1,n2,n3) = si_ifc(count,4);
                            asr_sum = asr_sum + si_ifc(count,4);
                            count = count + 1;
                        end
                    end
                end
            end
            % forces acoustic sum rule to "self" energy at ion [1,1,1]
            % (this makes w=0 at q=0 for the acoustic branches)
            phid(i,j,na,na,1,1,1) = phid(i,j,na,na,1,1,1) - asr_sum;
        end
    end
end

% more propoer MATLAB way to do it, but harder to impose acoustic sum rules:
% phid=reshape(si_ifc', 4, [], 36);
% phid=phid(4,2:end,:); 
% phid=reshape(phid,nr3,nr2,nr1,nat,nat,3,3);
% phid=permute(phid,[7:-1:1]);   

fclose(fp);
