%%
addpath(genpath('~/Dropbox/xdscode/'))
cd ~/Dropbox/xdscode/tests/

ry_to_thz           = 3289.842;
amu_to_ry           = 911.4442;
bohr_angstr         = 0.52917720859;
ry_to_ev            = 13.6057;

%%
coeffs = parse_scattfact_file('Si', '../scattering/scattfact/scattfact.txt')
as = coeffs(1:2:7);
bs = coeffs(2:2:8)/4.0;
c = coeffs(9);

%%
ff = make_form_factor('Ga')

%% 
xrayff = make_xrayscatt_factor('Ga',10500)

%%
elements = ['V','O'];
ii=1;
for el = elements
    g.form_factors{ii} = make_xrayscatt_factor(el, 9500)
    ii = ii+1;
end

%% 
s = load('vo2m1.mat','geometry');
% s = load('vo2rutile.mat','geometry');
geometry = s.geometry;

for ii = 1:length(geometry.element)
    el = geometry.element{ii}
    geometry.form_factors{ii} = make_xrayscatt_factor(el, 9500);
end

%%