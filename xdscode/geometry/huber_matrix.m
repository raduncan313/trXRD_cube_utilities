function Rot = huber_matrix(phi, theta, chi)
    yaxis = [0., -1., 0.]';
    zaxis = [0., 0., 1.]';
    M = rotationmat3D(theta, yaxis);
    xaxis = M*[1., 0., 0.]';
    Rot = rotationmat3D(chi, xaxis)*M*rotationmat3D(phi, zaxis);
end
