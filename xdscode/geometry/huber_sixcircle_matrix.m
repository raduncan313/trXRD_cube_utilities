function Rot = huber_sixcircle_matrix(phi, theta, chi, mu)
    yaxis = [0., 1., 0.]';
    zaxis = [0., 0., 1.]';
    xaxis = [1., 0., 0.]';

    % order of rotations: phi(y), chi(x), eta(y)
    
    Mtheta = rotationmat3D(theta, -yaxis);
    Mchi   = rotationmat3D(chi  , xaxis);
%     Mphi   = rotationmat3D(phi  , zaxis);
    Mphi   = rotationmat3D(phi  , -yaxis);
    Mmu = rotationmat3D(mu  , zaxis);
    Rot = Mmu*Mtheta*Mchi*Mphi;  % *rotationmat3D(-90, [1,0,0])  % <-- right multiply to compare with `huber_matrix` angles

    
end
