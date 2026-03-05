function C = stiffnessMatrixLocalToGlobal(C_local, R)
% Rotates elasticity tensor(s) from local to global coordinates
%
% Inputs
%   C_local : local stiffness matrix/matrices in Voigt notation
%             size = 6x6 or 6x6xNg
%   R       : rotation matrix/matrices (direction cosines), size = 3x3xNg
%
% Output
%   C       : rotated stiffness matrix/matrices in global coordinates
%             size = 6x6xNg
%
% Conventions
%   - Voigt order is assumed to be [11 22 33 23 13 12].
%   - Bond matrix below is consistent with this order.
%   - R must map local components to global components (local -> global).
%     If your R is global -> local, pass permute(R,[2 1 3]) instead.
%
% Notes
%   - Works for any symmetry class, as long as C_local uses the same Voigt convention.
%   - If C_local is 6x6, the same local tensor is used for all grains.
%   - If C_local is 6x6xNg, one local tensor per grain is used.
%
% IMPORTANT:
% R must map local crystal coordinates to global/lab coordinates.
% If Euler angles are interpreted in the opposite sense, use R = permute(R,[2 1 3]).

% Input checks
if ismatrix(R)
    if ~isequal(size(R), [3 3])
        error('R must be 3x3 or 3x3xNg.');
    end
    R = reshape(R, 3, 3, 1);
elseif ndims(R) == 3
    if size(R,1) ~= 3 || size(R,2) ~= 3
        error('R must be 3x3xNg.');
    end
else
    error('R must be 3x3 or 3x3xNg.');
end

% number of grains (can be 1)
Ng = size(R,3);

if ismatrix(C_local)
    if ~isequal(size(C_local), [6 6])
        error('C_local must be 6x6 or 6x6xNg.');
    end
    singleGrain = true;
elseif ndims(C_local) == 3
    if size(C_local,1) ~= 6 || size(C_local,2) ~= 6
        error('C_local must be 6x6 or 6x6xNg.');
    end
    if size(C_local,3) ~= Ng
        error('If C_local is 6x6xNg, its 3rd dimension must match size(R,3).');
    end
    singleGrain = false;
else
    error('C_local must be 6x6 or 6x6xNg.');
end


% Bond transformation matrix M (Voigt order [11 22 33 23 13 12])
M = zeros(6,6,Ng);

% Top-left block (normal <- normal)
M(1:3,1:3,:) = R(1:3,1:3,:).^2;

% Top-right block (normal <- shear)
M(1,4,:) = 2*R(1,2,:).*R(1,3,:);
M(1,5,:) = 2*R(1,3,:).*R(1,1,:);
M(1,6,:) = 2*R(1,1,:).*R(1,2,:);

M(2,4,:) = 2*R(2,2,:).*R(2,3,:);
M(2,5,:) = 2*R(2,3,:).*R(2,1,:);
M(2,6,:) = 2*R(2,1,:).*R(2,2,:);

M(3,4,:) = 2*R(3,2,:).*R(3,3,:);
M(3,5,:) = 2*R(3,3,:).*R(3,1,:);
M(3,6,:) = 2*R(3,1,:).*R(3,2,:);

% Bottom-left block (shear <- normal)
M(4,1,:) = R(2,1,:).*R(3,1,:);
M(4,2,:) = R(2,2,:).*R(3,2,:);
M(4,3,:) = R(2,3,:).*R(3,3,:);

M(5,1,:) = R(3,1,:).*R(1,1,:);
M(5,2,:) = R(3,2,:).*R(1,2,:);
M(5,3,:) = R(3,3,:).*R(1,3,:);

M(6,1,:) = R(1,1,:).*R(2,1,:);
M(6,2,:) = R(1,2,:).*R(2,2,:);
M(6,3,:) = R(1,3,:).*R(2,3,:);

% Bottom-right block (shear <- shear)
M(4,4,:) = R(2,2,:).*R(3,3,:) + R(2,3,:).*R(3,2,:);
M(4,5,:) = R(2,1,:).*R(3,3,:) + R(2,3,:).*R(3,1,:);
M(4,6,:) = R(2,2,:).*R(3,1,:) + R(2,1,:).*R(3,2,:);

M(5,4,:) = R(1,2,:).*R(3,3,:) + R(1,3,:).*R(3,2,:);
M(5,5,:) = R(1,3,:).*R(3,1,:) + R(1,1,:).*R(3,3,:);
M(5,6,:) = R(1,1,:).*R(3,2,:) + R(1,2,:).*R(3,1,:);

M(6,4,:) = R(1,2,:).*R(2,3,:) + R(1,3,:).*R(2,2,:);
M(6,5,:) = R(1,3,:).*R(2,1,:) + R(1,1,:).*R(2,3,:);
M(6,6,:) = R(1,1,:).*R(2,2,:) + R(1,2,:).*R(2,1,:);

Mt = permute(M, [2 1 3]); % transpose of M

% Note: pagemtimes is introduced after version 2022b
hasPagemtimes = ~isempty(which('pagemtimes'));

% Rotate stiffness tensor(s)
if hasPagemtimes
    % Works for both:
    %   C_local = 6x6      (same tensor for all grains)
    %   C_local = 6x6xNg   (one tensor per grain)
    C = pagemtimes(M, pagemtimes(C_local, Mt));
else
    C = zeros(6,6,Ng);
    if singleGrain
        for j = 1:Ng
            C(:,:,j) = M(:,:,j) * C_local * Mt(:,:,j);
        end
    else
        for j = 1:Ng
            C(:,:,j) = M(:,:,j) * C_local(:,:,j) * Mt(:,:,j);
        end
    end
end

end
