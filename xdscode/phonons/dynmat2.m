function DD2 = dynmat2(Qflat, ifc_struct, masses)
%%
twopii = -2*pi*1j;
ifc=ifc_struct.ifc;
multiplicities=ifc_struct.multiplicities;
alltaudiffs=ifc_struct.alltaudiffs;

m  = length(multiplicities);
nq = size(Qflat, 1);
nat= ifc_struct.num_at;
nuc= ifc_struct.num_uc;

exparg = exp(twopii*Qflat*alltaudiffs);
exparg2 = (zeros(nq, m));
i1=1;
for ii=1:m
    i2 = i1+multiplicities(ii)-1;
    exparg2(:,ii) = sum(exparg(:,i1:i2),2)/multiplicities(ii);
    i1 = i2+1;
end
exparg2 = reshape(exparg2, nq, nat, nat, nuc);
% exparg2 = reshape(exparg2', nuc, nat, nat, nq);

%% multiply force matrices by phases and normalize by sqrt(mass) then sum
DD2 = (zeros(3*nat,3*nat,nq));
for nb=1:nat
    for na=1:nat
        sifc=reshape(ifc(:,:,na,nb,:),3*3,nuc);
%         exparg4 = squeeze(exparg2(:,na,nb,:)).';
        exparg4 = reshape(exparg2(:,na,nb,:),[nq, nuc]).';
        sqrt_mamb = sqrt(masses(na)*masses(nb));
DD2( (3*(na-1)+1):3*na, (3*(nb-1)+1):3*nb, :) = reshape(sifc*exparg4,3,3,nq)/sqrt_mamb;
    end
end

end