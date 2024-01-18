function plotunitcell_withphonon(at, tau, vv, varargin)
    plotunitcell(at);
    hold on
    axis equal vis3d; 
    nat = size(tau,2);
    vv = reshape(vv, [3, nat]);
    
    scatter3(tau(1,:), tau(2,:), tau(3,:),75,'filled','b');
    
    if isempty(varargin)
        quiver3(tau(1,:), tau(2,:), tau(3,:), vv(1,:), vv(2,:), vv(3,:),'linewidth',1.5);
    else
        quiver3(tau(1,:), tau(2,:), tau(3,:), vv(1,:), vv(2,:), vv(3,:), varargin);
    end
end
