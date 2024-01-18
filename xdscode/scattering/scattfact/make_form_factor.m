function form_factor = make_form_factor(elem)
% coeffs = parse_scattfact_file(elem, '../scattering/scattfact/scattfact.txt');
coeffs = parse_scattfact_file(elem, [fileparts(which('make_form_factor')),'/scattfact.txt']);
as = coeffs(1:2:7);
bs = coeffs(2:2:8)/4.0;
c = coeffs(9);

form_factor.as = as;
form_factor.bs = bs;
form_factor.c = c;

end

