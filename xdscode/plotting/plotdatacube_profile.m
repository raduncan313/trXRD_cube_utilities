function plotdatacube_profile(varargin)
% function plotdatacube(t, data, geometry, caxis, box)
% 
thisf=gcf();
[t, data, geometry, caxis_, box] = parse_inputs(varargin{:});
if isempty(t)
    t = 1:size(data,3);
end

nt = length(t);
nticks = floor(nt/4);

imgax = axes('Position',[0.65 0.65 0.28 0.28]);
imobj = imagesc(imgax, data(:,:,1));
caxis(caxis_);

% plot time slice when using datatip
% fh = figure(); clf

% set the datacursor callback function
dcm_obj = datacursormode(thisf);
set(dcm_obj, 'DisplayStyle', 'window')
set(0,'CurrentFigure',thisf);

% create a toolbar
tbh = uitoolbar(thisf);
% Add a push tool to the toolbar
img2 = zeros(16,16,3);
tth = uipushtool(tbh,'CData',img2,'Separator','on',...
           'TooltipString','Select ROIs from these cube data',...
           'HandleVisibility','off');

nrois=1; % TODO: for now
tth.ClickedCallback = {@plotroi, t, data};

if ~isempty(geometry)
    geometry.imageNz = size(data,1);
    geometry.imageNy= size(data,2);
    m = geometry.imageNz;
    n = geometry.imageNy;
    [Q, QQ, allK_q, allhkl] = generate_reduced_wavevectors(geometry);
    plotbzcontours(geometry, Q, QQ, allhkl);
    bg = geometry.primvects';
    lambda0 = geometry.lambda0;
    theta   = geometry.theta;
    chi     = geometry.chi;
    phi0    = geometry.phi;
    SamRot 	= geometry.SamRot;
    Rot = huber_matrix(phi0, theta, chi)*SamRot;
%     MM = inv(bg)*inv(Rot);
    Q = reshape(Q, [m, n, 3]);
    QQ= reshape(QQ, [m, n, 3]);
    allhkl= reshape(allhkl, [m, n, 3]);
    set(dcm_obj,'UpdateFcn',{@myupdatefn_wgeom, ...
        {QQ, Q, allhkl, t, data, lambda0, bg, Rot}, thisf, box})
else
    set(dcm_obj,'UpdateFcn',{@myupdatefn_nogeom, t, data, thisf, box})
end

% imobj
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
set(hjSlider, 'StateChangedCallback', {@sliderCallback, imobj, t, data});

ii = get(hjSlider,'Value');
% geometry.phi = t(ii);

% NOTE: only pass `pos` then find QQ, hkl's only at `pos` duh!
1;

end %function


% %=====================================================================

function sums=selectroi(data)
    sums = zeros(size(data,3),1);
    bw = roipoly;
    for ii=1:size(data,3)
        sums(ii) = sum(sum(bw.*squeeze(data(:,:,ii))));
    end
end


% %=====================================================================
function td=plotroi(a, b, t, data)
%     sums=zeros([nrois, size(data,3)]);
%     for roict = 1:nrois
%         sums(roict,:) = selectroi(data);
%     end
    
    td = selectroi(data);
    td = td-min(td(:));
    fh=figure(3);
    subplot(2,1,1)
    title('ROI')
    plot(t,td)
    subplot(2,1,2)

    id = t > 0.;
    td = td(id);
    [ws, Xw] = tdsfft(t(id), td(:));
    xw = abs(Xw);
    xw = xw./max(xw);
    plot(ws, xw)
end

% TODO: needs updating...
% %=====================================================================
function plotbzcontours(geometry, Q, QQ, allhkl)
    Nz = geometry.imageNz;
    Ny = geometry.imageNy;

    Q=reshape(Q,[Nz, Ny, 3]);
    QQ=reshape(QQ,[Nz, Ny, 3]);

    %% FIXME: use gradient and scatter() to plot this
    hold on
    rr=abs(reshape(sqrt(2)*allhkl(:,1) - pi*allhkl(:,2) + exp(1)*allhkl(:,3), [Nz, Ny]));
    contour(rr,30,'color','w');
    allhkl=reshape(allhkl,[Nz, Ny, 3]);
    k0 = [1,0,0]'/geometry.lambda0;
    lambda0=geometry.lambda0;
    bg = geometry.primvects';

    theta   = geometry.theta;
    chi     = geometry.chi;
    phi0    = geometry.phi;
    SamRot 	= geometry.SamRot;
    Rot = huber_matrix(phi0, theta, chi)*SamRot;
    MM = inv(bg)*inv(Rot);
    xi = 0; yj = 0;

    %% use ginput to set BZ labels on graph
    while (xi < Ny & yj < Nz)
    %	figure(datafig);
        [x,y]=ginput(1);
        xi = round(x);	
        yj = round(y);
        try
            hkl = squeeze(allhkl(yj,xi,:));
            qq=squeeze(QQ(yj,xi,:));
            q=squeeze(Q(yj,xi,:));
            kp = k0+q;
            mu = atan2d(kp(2), kp(1));
            del = 90-acosd(kp(3)*geometry.lambda0);
            text(x,y,sprintf('(%i %i %i)',hkl),'FontWeight','bold','Color','w','fontsize',14);

            disp(['(i, j) = (', num2str([xi, yj]),')' ])
            disp([sprintf('(%i %i %i)\t',hkl), ...
                'ThetaB = ', num2str(asind(lambda0*norm(q)/2)),...
                '    mu = ', num2str( mu), ...
                '    del = ', num2str(del) ])
            disp(['q = [',num2str((MM*qq)'), ']'])
            disp(['Q = [',num2str((MM*q)'), ']'])

        catch
           disp('Done.'); 
        end
        %}
    end
    %}
    1;

end

% %=====================================================================
function sliderCallback(hObject, hEventData, imobj, t, img)
    ii = get(hObject,'Value');
    set(imobj, 'cdata', img(:,:,ii))
    title(['t = ', num2str(t(ii)), '  i = ', num2str(ii) ]);
end

% %=====================================================================
function txt = myupdatefn_wgeom(~,event_obj,arrbundle, origf, box)  %
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
% origf = gcf();
fh = figure(2);
% I = get(event_obj, 'DataIndex');
q = squeeze(Q(pos(2), pos(1), :));
qq= squeeze(QQ(pos(2), pos(1), :));
hkl=squeeze(allhkl(pos(2), pos(1), :));
% MM = inv(bg);
k0 = [1,0,0]'/lambda0;
kp = k0+q(:);
mu = atan2d(kp(2), kp(1));
del = 90-acosd(kp(3)*lambda0);
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       sprintf('hkl: %i, %i, %i',hkl),...
       sprintf('mu: %g    del: %g', mu, del),...      %['mu: ', num2str(mu), '\tdel: ', num2str(del)]
       sprintf('ThetaB: %g    inv(d): %g', asind(lambda0*norm(q)/2), 1/norm(q)),...
       sprintf('q : %1.3g  %1.3g  %1.3g', bg\qq),...
       sprintf('Q : %1.3g  %1.3g  %1.3g', (Rot*bg)\q)};
       
td=get_time_trace_from_cube(dd, pos, box);
% figure(fh)
set(0,'CurrentFigure',fh);
subplot(2,1,1)
plot(t,td)
title(['x, y = ', num2str(pos(1)), ', ' num2str(pos(2))])
subplot(2,1,2)
[ws, Xw] = tdsfft(t, td);
ids = ws>0.2 & ws < 10;
xw = abs(Xw);
% xw = xw-mean(xw(end-20:end));
xw = xw./max(xw(ids));
plot(ws, xw)
xlim([0, 10])
ylim([0, 1.2])
% txt = {['X: ',num2str(pos(1))],...
%        ['Y: ',num2str(pos(2))]};
set(0,'CurrentFigure',origf);

end

%=====================================================================
function txt = myupdatefn_nogeom(~,event_obj,t,dd, origf, box)  %
pos = get(event_obj,'Position');
% origf = gcf();  % FIXME: could find from event_obj?
td=get_time_trace_from_cube(dd, pos, box);
fh=figure(2);
set(0,'CurrentFigure',fh);
subplot(2,1,1)
title(['x, y = ', num2str(pos(1)), ', ' num2str(pos(2))])
plot(t,td,'.-')
subplot(2,1,2)

id = t > 0.;
% [ws, Xw] = tdsfft(t, td');
td = td(id);
[ws, Xw] = tdsfft(t(id), td(:));
% ids = ws>0.2 & ws < 10;
xw = abs(Xw);
% xw = xw-mean(xw(end-20:end));
xw = xw./max(xw);
plot(ws, xw)
% xlim([0, 10])
% ylim([0, 1.2])

% get the intensity data
cdat=event_obj.Target.CData(pos(2),pos(1));

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
      sprintf('CData : %2.2E', cdat) };

set(0,'CurrentFigure',origf);
end


% %=====================================================================
function [t, data, geometry, caxis_, box] = parse_inputs(varargin)

narginchk(2,5)

box = [0,0];
caxis_ = [-Inf, Inf];
geometry = [];

switch nargin
    case 2,                     % plotdatacube(t, data)
        t = varargin{1};
        data = varargin{2};
    case 3,                     % plotdatacube(t, data, geometry)
        t = varargin{1};
        data = varargin{2};
        geometry = varargin{3};
    case 4,                     % plotdatacube(t, data, geometry, box)
        t = varargin{1};
        data = varargin{2};
        geometry = varargin{3};
        caxis_ = varargin{4};
        
    case 5,                     % plotdatacube(t, data, geometry, box, caxis)
        t = varargin{1};
        data = varargin{2};
        geometry = varargin{3};
        caxis_ = varargin{4};
        box = varargin{5};

    otherwise,
        error(message('plotdatacube:invalidInputs'));

end % switch
end 
% 
