function [intens, ww, vv] = tds_intensity_onq(Qflat, geometry, ifc)

ry_to_thz=3289.842;
amu_to_ry=911.4442;

masses = geometry.masses*amu_to_ry;

DD2 = dynmat2(Qflat, ifc, masses);
[ww, vv] = phonon_freqs_vecs(DD2);
intens=tdsI1(geometry, Qflat, ww, vv);
ww = ww*ry_to_thz;
end
