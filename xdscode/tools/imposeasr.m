function ifc = imposeasr(ifc)
	nat= size(ifc.ifc, 3);
	nuc= size(ifc.ifc, 5);
	for x=1:3
		for y=1:3
            for na = 1:nat
        		asrsum=0;
        		for nb = 1:nat
	        		for uc = 1:nuc
	        			asrsum = asrsum + ifc.ifc(x,y,nb,na,uc);
	        		end
	        	end
	        	ifc.ifc(x,y,na,na,1) = ifc.ifc(x,y,na,na,1) - asrsum;
        	end
		end
    end
end
