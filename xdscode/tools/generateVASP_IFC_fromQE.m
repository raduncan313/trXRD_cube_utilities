%% convert IFC from QE to VASP shape
% [at, a0, nr1, nr2, nr3, phid] = parse_ifc_file('ge444.fc');
% the FCC unit cell in QE is defined as this:
%       v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)
[at, tau, masses, nr1, nr2, nr3, phid] = parse_ifc_file('si444.fc');
% [at, tau, masses, nr1, nr2, nr3, phid] = parse_ifc_file('ge444.fc');
%% assumes 4x4x4 supercell
A_s =   at*[nr1 0 0
            0 nr2 0
            0 0 nr3 ];

%% generate list of neighbors on an nr1 x nr2 x nr3  supercell (tau_s)
nat = size(tau,2);
nuc = nr1*nr2*nr3;
sphi=size(phid);
natsc = (prod(sphi(4:end)));
tau_s = zeros(3,natsc);
count = 0;
for na=1:nat
    for n3=1:nr3
        for n2=1:nr2
            for n1=1:nr1
%                 r = at*[(n1-1) (n2-1) (n3-1)]' + tau(:,na);
                r = at*[n1-1 n2-1 n3-1]' + tau(:,na);
                count = count + 1;
                tau_s(:,count) = r;
            end
        end
    end
end

%% test generate tau_diffs
ct = 1;
tau_diffs = [];
% for n2=1:nuc
n2 = 1; % <---- how to choose this?
	for n1=1:nuc
        for nb=1:nat
			for na=1:nat
				RR = tau_s(:,n1-1+nuc*(na-1)+1)-tau_s(:,n2-1+nuc*(nb-1)+1);
				shortvecs=findinwscell(RR, A_s, 1e-12);
                tau_diffs{ct} = A_s*shortvecs;
                ct = ct+1;
			end
        end
	end
% end

%% check: visualize the tau_diffs genearated from tau_s, should all be in the WS cell
figure(5); clf
plotunitcell(A_s); 
scatter3(tau_s(1,:), tau_s(2,:), tau_s(3,:),'b')
box on
plotBZ(A_s);
aa=cell2mat(tau_diffs);
scatter3(aa(1,:), aa(2,:), aa(3,:),'r');
axis tight vis3d

%% the MEAT: generate VASP forces and corresponding tau_diffs for dynmat2():
%% generate both tau_diff and VASP IFC
nat = size(tau,2);
nuc = nr1*nr2*nr3;
FC = zeros(3,3,nat*nat*nuc);
tau_diffs=[];
ct = 1;
for n3=1:nr3
    for n2=1:nr2
        for n1=1:nr1
            for nb=1:nat
                for na=1:nat
                    r = at*[n1-1 n2-1 n3-1]' + (tau(:,na)-tau(:,nb));
                    shortvecs=findinwscell(r, A_s, 1e-6);
                    tau_diffs{ct} = A_s*shortvecs;
                    FC(:,:,ct) = phid(:,:,na,nb,n1,n2,n3);
                    ct = ct+1;
                end
            end
        end
    end
end
FC = reshape(FC, 3,3, nat, nat, nuc);

%%
ifc.ifc = FC;
ifc.tau_diffs=tau_diffs;  % this is not really used later
for k=1:length(tau_diffs)
    multiplicities(k) = size(tau_diffs{k},2);
end
ifc.multiplicities=multiplicities;
ifc.alltaudiffs=cell2mat(tau_diffs); % this is the one we use, flattened
ifc.num_at = nat;
ifc.num_uc = nuc;

%% save geometry and IFC file for use with dynmat2():
load('si.mat','geometry');
% load('ge.mat','geometry');

%%
save('si.mat','geometry','ifc');
disp('Saved `geometry` and `ifc` to si.mat')

