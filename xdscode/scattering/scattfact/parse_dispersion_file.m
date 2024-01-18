function [f1, f2] = parse_dispersion_file(elem, energy)
	file = [fileparts(which('parse_dispersion_file')), '/sf/' lower(elem), '.nff'];
	fid = fopen(file);
	hdr = fgetl(fid); % discard header
	str = fgetl(fid);
% 	en_prev, f1_prev, f2_prev = map(x->parse(Float64, x), split(str))
    tmp = str2double(split(str));
    en_prev = tmp(1); 
    f1_prev = tmp(2);
    f2_prev = tmp(3);
	while ~feof(fid)
		str = fgetl(fid);
% 		en_next, f1_next, f2_next = map(x->parse(Float64, x), split(str))
        tmp = str2double(split(str));
        en_next = tmp(1); 
        f1_next = tmp(2);
        f2_next = tmp(3);
		if en_prev <= energy && en_next > energy
			den = (en_next - en_prev);
			df1 = f1_next - f1_prev ;
			df2 = f2_next - f2_prev;
			f1 = f1_prev + df1/den*(energy-en_prev);
			f2 = f2_prev + df2/den*(energy-en_prev);
			return;
		end
        en_prev = en_next;
        f1_prev = f1_next;
        f2_prev = f2_next;
	end
	warning("Scattfactors:Energy out of range!")
end
