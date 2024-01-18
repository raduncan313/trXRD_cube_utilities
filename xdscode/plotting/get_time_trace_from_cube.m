function td = get_time_trace_from_cube(dd, pos, box)
% gets a 1D trace from the cube data `dd` at coordinates given by `pos`
% if `box` is not empty gives an average over a box size 2box(1)+1 x 2box(2)+1 
%
j= pos(1);
i= pos(2);
n = size(dd,1);
m = size(dd,2);
if isempty(box)
    td = squeeze(dd(i, j, :));
else
    i1 = i-box(1);
    i2 = i+box(1);
    i1 = max(1, i1);
    i2 = min(i2, n);

    j1 = j-box(2);
    j2 = j+box(2);
    j1 = max(1, j1);
    j2 = min(j2, m);
    td = squeeze(mean(mean(dd(i1:i2, j1:j2, :),1),2));
end
