function ffxray = xray_form_factors(geometry, Qflat)
    nq  = size(Qflat,1);
%     nat = geometry.num_atoms;
    nat = length(geometry.form_factors);  % should be == geometry.num_atoms;
    q2 = sum(Qflat.^2,2);
    ffxray = zeros(nq, nat);
    for na=1:nat
        fxray = geometry.form_factors{na};
        fatom = fxray.ff;
        Z  = sum(fatom.as)+fatom.c; % this is the atomic number
        ffxray(:, na) = atomic_form_factor(q2, fatom) + fxray.f1 + 1j*fxray.f2 - Z;
    end
end

function fout = atomic_form_factor(q2, fin)
%
	fout = fin.c;
    for i=1:length(fin.as)
		fout = fout + fin.as(i)*exp(-fin.bs(i).*q2);
    end
end
