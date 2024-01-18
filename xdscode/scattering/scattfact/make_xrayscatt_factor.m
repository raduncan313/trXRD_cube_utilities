function xray_factor = make_xrayscatt_factor(elem,energy)
    ff = make_form_factor(elem);
    [f1, f2] = parse_dispersion_file(elem, energy);
    xray_factor.ff = ff;
    xray_factor.energy = energy;
    xray_factor.f1 = f1;
    xray_factor.f2 = f2;

end

