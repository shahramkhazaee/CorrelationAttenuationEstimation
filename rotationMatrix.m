function R = rotationMatrix(phi1, Phi, phi2)
% Bunge style ZXZ rotation, implemented as Rz(phi1)*Rx(Phi)*Rz(phi2)
%
% Inputs
%   phi1, Phi, phi2 : scalars or vectors (Ngx1 or 1xNg). Angles in radians.
%                     Scalars are expanded to match the vector length.
%
% Output
%   R : 3x3xNg array of rotation matrices

    % force to column vectors
    phi1 = phi1(:);
    Phi  = Phi(:);
    phi2 = phi2(:);

    Ng = max([numel(phi1), numel(Phi), numel(phi2)]);

    % expand scalars
    if numel(phi1) == 1, phi1 = repmat(phi1, Ng, 1); end
    if numel(Phi)  == 1, Phi  = repmat(Phi,  Ng, 1); end
    if numel(phi2) == 1, phi2 = repmat(phi2, Ng, 1); end

    % final size check
    if ~(numel(phi1)==Ng && numel(Phi)==Ng && numel(phi2)==Ng)
        error('phi1, Phi, phi2 must be scalars or vectors of the same length.');
    end

    % reshape for page wise operations: 1x1xNg
    phi1 = reshape(phi1, 1, 1, Ng);
    Phi  = reshape(Phi,  1, 1, Ng);
    phi2 = reshape(phi2, 1, 1, Ng);

    z = zeros(size(phi1));
    o = ones(size(phi1));

    R1 = [cos(phi1) -sin(phi1) z;
          sin(phi1)  cos(phi1) z;
          z          z         o];

    z = zeros(size(Phi));
    o = ones(size(Phi));

    R2 = [o  z         z;
          z  cos(Phi) -sin(Phi);
          z  sin(Phi)  cos(Phi)];

    z = zeros(size(phi2));
    o = ones(size(phi2));

    R3 = [cos(phi2) -sin(phi2) z;
          sin(phi2)  cos(phi2) z;
          z          z         o];

    % Note: pagemtimes is introduced after version 2022b
    hasPagemtimes = ~isempty(which('pagemtimes'));

    % Multiply page wise
    if hasPagemtimes
        R = pagemtimes(pagemtimes(R1, R2), R3);
    else
        % Fallback loop for older MATLAB
        R = zeros(3,3,Ng);
        for k = 1:Ng
            R(:,:,k) = R1(:,:,k) * R2(:,:,k) * R3(:,:,k);
        end
    end
    
end

