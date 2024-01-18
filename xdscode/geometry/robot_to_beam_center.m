function [xc,yc,dd] = robot_to_beam_center(x,y,z)
% robot for r138 (InGaAs with rings)
pixel_size = 0.1;
% x0 = 55.502;
% y0 = 159.591;
% z0 = 510.049;

x0 = 0.004;
y0 = 94.133;
z0 = 83.155;

det_dist = 117;
beam_center = [812, 1783];

dx = x-x0;
dy = y-y0;
dz = z-z0;
xc = dx/pixel_size + beam_center(1); % these are in pixel units
yc = dy/pixel_size + beam_center(2);
dd = dz + det_dist;

end
