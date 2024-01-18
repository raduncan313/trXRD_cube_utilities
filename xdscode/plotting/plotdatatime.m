function plotdatatime(t, data, varargin)
thisf=gcf();
box = [0,0];
rng = [-inf, inf];

if isempty(t)
    t = 1:size(data,3);
end
if nargin <= 4
    rng = varargin{1};
end
if nargin == 4
    box = varargin{2};
end

nt = length(t);
nticks = floor(nt/4);

imobj = imagesc(data(:,:,1), rng);
axis image
set(gca, 'xlimmode','manual',...
'ylimmode','manual',...
'zlimmode','manual',...
'climmode','manual',...
'alimmode','manual');

% colorbar

% setup java slider control
[jSlider, hSlider] = javacomponent('javax.swing.JSlider','South');
set(jSlider, 'Value',1, 'MajorTickSpacing', nticks, 'PaintLabels',true);
jSlider.setMinimum(1)
jSlider.setMaximum(nt)
hjSlider = handle(jSlider, 'CallbackProperties');
set(hjSlider, 'StateChangedCallback', {@sliderCallback, imobj, t, data,rng});

% plot time slice when using datatip
fh = figure(); clf

% set the datacursor callback function
dcm_obj = datacursormode(thisf);
set(dcm_obj, 'DisplayStyle', 'window')
set(dcm_obj,'UpdateFcn',{@myupdatefcn, t, data, fh, box})
1;
end


%
% %=====================================================================
% function params = parseInputs(varargin)
% 
% end
% 


function sliderCallback(hObject, hEventData, imobj, t, img, rng)
%     disp(hEventData)
    ii = get(hObject,'Value');
%     imagesc(img(:,:,ii), rng)
    set(imobj, 'cdata', img(:,:,ii))
    axis image
    title(['t = ', num2str(t(ii)), '  i = ', num2str(ii) ]);
    
end

%=====================================================================
function txt = myupdatefcn(~,event_obj,t,dd, fh, box)  %
pos = get(event_obj,'Position');
j= pos(1);
i= pos(2);
% n = size(dd,1);
% m = size(dd,2);
% if isempty(box)
%     td = squeeze(dd(i, j, :));
% else
%     i1 = i-box(1);
%     i2 = i+box(1);
%     i1 = max(1, i1);
%     i2 = min(i2, n);
% 
%     j1 = j-box(2);
%     j2 = j+box(2);
%     j1 = max(1, j1);
%     j2 = min(j2, m);
%     td = squeeze(mean(mean(dd(i1:i2, j1:j2, :),1),2));
% end
td=get_time_trace_from_cube(dd, pos, box);
figure(fh)
subplot(2,1,1)
plot(t,td)
subplot(2,1,2)
[ws, Xw] = tdsfft(t, td');
ids = ws>0.2 & ws < 10;
xw = abs(Xw);
xw = xw-mean(xw(end-20:end));
xw = xw./max(xw(ids));
plot(ws, xw)
xlim([0, 10])
ylim([0, 1.2])
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))]};
end