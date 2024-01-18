function [a, Rt, attypes_p, atnums_p]=readPOSCAR(file)

%% Read in Unit Cell Positions

fp=fopen(file,'r'); % open the specified POSCAR file

%for j=1:2
test = fgetl(fp);%read in and discard header
%disp(test)
%end

% read in this scaling factor they have for the lattice vectors
sf=fscanf(fp,'%f',[1,1]); 

test = fgetl(fp);%read in and discard header
%disp(test)
%
a=fscanf(fp,'%f',[3,3]); % read in the lattice vectors
fgetl(fp);

% scale the lattice vectors by the scaling factor
a = a*sf;

test =fgetl(fp);
attypes_p = regexp(test,'[a-zA-Z]+','match');
% unqionms = strsplit(fgetl(fp), ' '); % read in names of elements
% unqionms(strcmp('',unqionms))=[];

%
atnums_p= fscanf(fp,'%f',[length(attypes_p),1]); % number of each element

Nat = sum(atnums_p); % number of ions in unit cell

for j = 1:2
    fgetl(fp); % disregard this line
    %disp(test)
end

% read in ion positions
% these are ion positions in terms of the lattice vectors 
Rt=fscanf(fp,'%f',[3,sum(atnums_p)]);
%}
fclose(fp); % close file

% now the 'translated' ion positions, which are in terms of x y and z,
% units are angstroms?
% Rt = a*R;

end
