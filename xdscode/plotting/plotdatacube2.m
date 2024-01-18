function plotdatacube2(varargin)
thisf=gcf();

[t, data, geometry, caxis_, box] = parse_inputs(varargin{:});

nt = length(t);
nticks = floor(nt/4);
imobj = imagesc(data(:,:,1));
axis image
set(gca, 'xlimmode','manual',...
'ylimmode','manual',...
'zlimmode','manual',...
'climmode','manual',...
'alimmode','manual');
caxis(caxis_);

% setup java slider control
[jSlider, hSlider] = javacomponent('javax.swing.JSlider','South');
set(jSlider, 'Value',1, 'MajorTickSpacing', nticks, 'PaintLabels',true);
jSlider.setMinimum(1)
jSlider.setMaximum(nt)
hjSlider = handle(jSlider, 'CallbackProperties');
% set(hjSlider, 'StateChangedCallback', {@sliderCallback, imobj, t, data,caxis_});
if isempty(geometry)
    set(hjSlider, 'StateChangedCallback', {@sliderCallback_nogeometry, imobj, t, data});
else
%     geometry.imageNy = geometry.imageNy/2;
%     geometry.imageNz = geometry.imageNz/2;
    set(hjSlider, 'StateChangedCallback', {@sliderCallback_wgeometry, imobj, t, data, geometry});
    
%     allhkl = hjSlider(sliderCallback_wgeometry(imobj, t, data, geometry));
end
ii = get(hjSlider,'Value');
% geometry.phi = t(ii);

% plot time slice when using datatip
fh = figure(2); clf

% set the datacursor callback function
dcm_obj = datacursormode(thisf);
set(dcm_obj, 'DisplayStyle', 'window')
% set(dcm_obj,'UpdateFcn',{@myupdatefcn, t, data, fh, box})
% set(dcm_obj,'UpdateFcn',{@myupdatefcn2, {QQ, Q, allhkl, t, data, lambda0, bg, Rot}, fh, box})
set(dcm_obj,'UpdateFcn',{@myupdatefcn3, geometry, t, data, fh, box})

% NOTE: only pass `pos` then find QQ, hkl's only at `pos` duh!
1;

end %function


%
% %=====================================================================
function [t, data, geometry, caxis_, box] = parse_inputs(varargin)

narginchk(2,5)

box = [1,1];
caxis_ = [-Inf, Inf];
geometry = [];

switch nargin
    case 2,                     % plotdatacube(t, data)
        t = varargin{1};
        data = varargin{2};
        if isempty(varargin{1})
            t = 1:size(data,3);
        end
    case 3,                     % plotdatacube(t, data, geometry)
        t = varargin{1};
        data = varargin{2};
        if isempty(varargin{1})
            t = 1:size(data,3);
        end
        geometry = varargin{3};
    case 4,                     % plotdatacube(t, data, geometry, box)
        t = varargin{1};
        data = varargin{2};
        if isempty(varargin{1})
            t = 1:size(data,3);
        end
        geometry = varargin{3};
        box = varargin{4};
        
    case 5,                     % plotdatacube(t, data, geometry, box, caxis)
        t = varargin{1};
        data = varargin{2};
        if isempty(varargin{1})
            t = 1:size(data,3);
        end
        geometry = varargin{3};
        box = varargin{4};
        caxis_ = varargin{5};

    otherwise,
        error(message('plotdatacube:invalidInputs'));

end % switch
end % function parse_input
% 

% %=====================================================================
function sliderCallback(hObject, hEventData, imobj, t, img, rng)
%     disp(hEventData)
    ii = get(hObject,'Value');
%     imagesc(img(:,:,ii), rng)
    set(imobj, 'cdata', img(:,:,ii))
    axis image
    title(['t = ', num2str(t(ii)), '  i = ', num2str(ii) ]);
    
end

% %=====================================================================
function sliderCallback_wgeometry(hObject, hEventData, imobj, phis, img, geometry)
    ii = get(hObject,'Value');
    nz = size(img,1);
    ny = size(img,2);
    geometry.imageNz = floor(nz/4);
    geometry.imageNy = floor(ny/4);
    geometry.phi = phis(ii);
    bz = makebzmesh(geometry, [nz, ny]);
    set(imobj, 'cdata', img(:,:,ii).*(1+5*bz))
    title(['t = ', num2str(phis(ii)), '  i = ', num2str(ii) ]);
end

% %=====================================================================
function sliderCallback_nogeometry(hObject, hEventData, imobj, t, img)
    ii = get(hObject,'Value');
    set(imobj, 'cdata', img(:,:,ii))
    title(['t = ', num2str(t(ii)), '  i = ', num2str(ii) ]);
end

% %=====================================================================
function txt = myupdatefcn2(~,event_obj,arrbundle, fh, box)  %
QQ = arrbundle{1};
Q = arrbundle{2};
allhkl = arrbundle{3};
t = arrbundle{4};
dd = arrbundle{5};
lambda0 = arrbundle{6};
bg = arrbundle{7};
Rot = arrbundle{8};

% Customizes text of data tips
pos = get(event_obj,'Position');
% I = get(event_obj, 'DataIndex');
q = squeeze(Q(pos(2), pos(1), :))
qq= squeeze(QQ(pos(2), pos(1), :));
hkl=squeeze(allhkl(pos(2), pos(1), :));
MM = inv(bg)*inv(Rot);
k0 = [1,0,0]'/lambda0;
kp = k0+q(:);
mu = atan2d(kp(2), kp(1));
del = 90-acosd(kp(3)*lambda0);
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       sprintf('hkl: %i, %i, %i',hkl),...
       sprintf('mu: %g    del: %g', mu, del),...  %['mu: ', num2str(mu), '\tdel: ', num2str(del)]
       sprintf('ThetaB: %g    inv(d): %g', asind(lambda0*norm(q)/2), 1/norm(q)),...
       sprintf('q : %1.3g  %1.3g  %1.3g', MM*qq),...
       sprintf('Q : %1.3g  %1.3g  %1.3g', inv(bg)*q)};
       
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
% txt = {['X: ',num2str(pos(1))],...
%        ['Y: ',num2str(pos(2))]};

end

%=====================================================================
function txt = myupdatefcn(~,event_obj,t,dd, fh, box)  %
pos = get(event_obj,'Position');
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


%=====================================================================
function txt = myupdatefcn3(~,event_obj,geometry, t, imgdata, fh, box)  %
pos = get(event_obj,'Position');
geometry.imageNz = size(imgdata,1);
geometry.imageNy = size(imgdata,2);
bg = geometry.primvects';
theta   = geometry.theta;
chi     = geometry.chi;
phi0    = geometry.phi;
SamRot 	= geometry.SamRot;
Rot = huber_matrix(phi0, theta, chi)*SamRot;
lambda0 = geometry.lambda0;

Q=det_kspace_proj(geometry);
q = squeeze(Q(pos(2), pos(1), :));
hkl = findclosesthkl2(bg, (Rot\q)')'
qq= Rot\q - bg*hkl

MM = inv(bg)*inv(Rot);
k0 = [1,0,0]'/lambda0;
kp = k0+q(:);
mu = atan2d(kp(2), kp(1));
del = 90-acosd(kp(3)*lambda0);
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       sprintf('hkl: %i, %i, %i',hkl),...
       sprintf('mu: %g    del: %g', mu, del),...  %['mu: ', num2str(mu), '\tdel: ', num2str(del)]
       sprintf('ThetaB: %g    inv(d): %g', asind(lambda0*norm(q)/2), 1/norm(q)),...
       sprintf('q : %1.3g  %1.3g  %1.3g', MM*qq),...
       sprintf('Q : %1.3g  %1.3g  %1.3g', inv(bg)*q)};

td=get_time_trace_from_cube(imgdata, pos, box);
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
% txt = {['X: ',num2str(pos(1))],...
%        ['Y: ',num2str(pos(2))]};
end
