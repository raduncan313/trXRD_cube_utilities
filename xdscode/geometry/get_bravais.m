function at=get_bravais(a, b, c, al, be, ga) 
	at = [
		a 0 0;
		b*cosd(ga) b*sind(ga) 0;
		c*cosd(be) c*(cosd(al)-cosd(be)*cosd(ga))/sind(ga) c*sqrt( 1 + 2*cosd(al)*cosd(be)*cosd(ga) - cosd(al)^2-cosd(be)^2-cosd(ga)^2 )/sind(ga)
		]';
end
