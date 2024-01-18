function [Sa, SR, atnums_s] = readSPOSCAR(file)

%% Read in Supercell Positions

fp=fopen(file,'r'); % open the specified POSCAR file

%for j=1:2 
test = fgetl(fp);%read in and discard header
%disp(test)
%end

% read in this scaling factor they have for the lattice vectors
Ssf=fscanf(fp,'%f',[1,1]); 

% test = fgetl(fp);%read in and discard header
%disp(test)
%

Sa=fscanf(fp,'%f',[3,3]);

Sa = Sa*Ssf;

test=fgetl(fp)
test=fgetl(fp)

% for j = 1:3
% ?    test=fgetl(fp);
    %disp(test)
% end

% atnums_s = fscanf(fp,'%f',[1,length(split(test))]);
atnums_s = str2double(regexp(test,'\d+','match'));

% test=fgetl(fp)
test=fgetl(fp) % Direct?

SR=fscanf(fp,'%f',[3,inf]);


fclose(fp);

% SRt = Sa*SR; % supercell ion positions in x y z coordinates

SNat = sum(atnums_s); % number of ions in supercell, must be == size(SR,2)
if SNat ~= size(SR, 2)
    warning('Inconsistent atomic numbers in SPOSCAR')
end

%SNuc = SNat/Nat; % number of unit cells
