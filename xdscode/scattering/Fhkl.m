function [allFhkl, thetaB, dis] = Fhkl(geometry, allhkl)
    if iscolumn(allhkl)
        allhkl = allhkl';
    end
	bg = geometry.primvects';
	at = geometry.realvecs';
	tau = at\geometry.basis';    % convert to primitive coords.

	fs = geometry.form_factors;
	lambda0 = geometry.lambda0;
	q2 = sum((bg*allhkl').^2, 1)';
    dis = 1/sqrt(q2);
	allFhkl=zeros(size(q2,1),1);

    for s=1:size(tau, 2)
        xff = xray_scattfactor(fs{s}, q2);
		allFhkl = allFhkl + xff(:).*exp(2*pi*1j*allhkl*tau(:,s)); 
	end
	thetaB = asind(lambda0*sqrt(q2)/2);
end

function s=formfactor(f, q2)
	s = f.c;
	for i=1:length(f.as)
		s = s + f.as(i)*exp(-f.bs(i).*q2);
    end
end

function xrayff=xray_scattfactor(f, q2)
    ffact = formfactor(f.ff, q2);
    Z = formfactor(f.ff, 0.0);
    xrayff = ffact + f.f1 + 1j*f.f2 - Z;
end
