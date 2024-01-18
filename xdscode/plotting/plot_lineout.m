function qs = plot_lineout(geometry,data,varargin)
bg = geometry.primvects';
[data1d, qs] = qslineout(geometry,data,varargin{1});

figure();
subplot(2,1,1)
plot(data1d)
xlabel('Path index')
ylabel('intensity')
subplot(2,1,2)
plot(qs*inv(bg)')
xlabel('Path index')
ylabel('hkl (rlu)')

end

