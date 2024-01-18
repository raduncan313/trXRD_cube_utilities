function coeffs=parse_scattfact_file(elem, file)
	fid = fopen(file);
    coeffs=[];
    str = fgetl(fid);
	while ~feof(fid)
		str = fgetl(fid);
		splitstr = split(str);
		if strcmp(splitstr{1}, elem)
            coeffs = str2double(splitstr(3:11))';
            fclose(fid);
            return
% 			return [parse(Float64, x) for x in splitstr[3:11] ]
            
		end
	end
    fclose(fid);
	warning("Element not found!")
    
end
