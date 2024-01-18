
function I_q=tdsI1(geometry, Qflat, ww, vv)
%%
    amu_to_ry=911.4442;
    ry_to_thz = 3289.8;
    K_to_ry = 6.333630685659031e-6;
    masses = geometry.masses*amu_to_ry;
	Ndim = size(ww,1);
    nat = Ndim/3;
	fxray = xray_form_factors(geometry, Qflat);
	DW_factors = geometry.DW_factors;
	T=geometry.temp*K_to_ry;
    I_q = zeros(size(ww));
	q2 = sum(Qflat.^2, 2);
    nq = size(Qflat,1);
    fexpfactor = fxray.*exp(-q2*DW_factors');
    fexpfactor = bsxfun(@rdivide, fexpfactor, sqrt(masses)');
    vv2=permute(vv, [2,1,3]);
    % repeated arrays for matching the final array shapes
    Qflat2=repmat(Qflat', nat,1);
    fexpfactor2=fexpfactor(:,ceil((1:Ndim)/3))';
    %%
    for nj = 1:Ndim
        I_q(nj,:) = sum(fexpfactor2.*Qflat2.*squeeze(vv2(nj,:,:,:)),1);
%         I_q(nj,:) = sum(fexpfactor2.*Qflat2.*reshape(vv2(nj,:,:,:),[],nq),1);
    end
    I_q = (1./ww)*1./(tanh(0.5*ww./T)).*abs(I_q).^2;

end
