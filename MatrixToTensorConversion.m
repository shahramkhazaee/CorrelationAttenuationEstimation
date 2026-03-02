function T = MatrixToTensorConversion(C, type)
% Converts a 6x6 elasticity matrix to a 3x3x3x3 tensor
%
%   T = Mat6toTens4(C)
%   T = Mat6toTens4(C, type)
%
% Inputs
%   C    : 6x6 matrix (elasticity in Voigt-like notation)
%   type : optional string
%          - omitted / '' : direct mapping using index order [11 22 33 23 13 12]
%          - 'Kelvin'     : applies sqrt(2) scaling on shear terms 
%                           (Kelvin/Mandel-like normalization)
%
% Output
%   T : 3x3x3x3 tensor
%
% Notes
%   - Index order is assumed to be [11, 22, 33, 23, 13, 12].
%   - This function assumes minor symmetries and maps both (i,j) and (j,i)
%     to the same Voigt index.

% Input checks
if ~ismatrix(C) || any(size(C) ~= [6 6])
    error('C must be a 6x6 matrix.');
end

if nargin < 2
    mode = 'direct';
else
    if ~(ischar(type) || isstring(type))
        error('type must be a character vector or string.');
    end
    type = char(string(type));
    if strcmpi(type, 'Kelvin')
        mode = 'scaled';
    else
        error('Unknown type "%s". Use ''Voigt'' or omit the argument.', type);
    end
end

% -----------------------------
% Pair -> Voigt index map
% Order: [11,22,33,23,13,12]
% -----------------------------
vmap = [1 6 5;
        6 2 4;
        5 4 3];

T = zeros(3,3,3,3);

for i = 1:3
    for j = 1:3
        I = vmap(i,j);
        for k = 1:3
            for l = 1:3
                J = vmap(k,l);
                if strcmp(mode, 'direct')
                    T(i,j,k,l) = C(I,J);
                else
                    % Kelvin/Mandel-like normalized shear handling
                    if I <= 3 && J <= 3
                        T(i,j,k,l) = C(I,J);
                    elseif (I <= 3 && J >= 4) || (I >= 4 && J <= 3)
                        T(i,j,k,l) = C(I,J) / sqrt(2);
                    else % I>=4 && J>=4
                        T(i,j,k,l) = C(I,J) / 2;
                    end
                end
            end
        end
    end
end

end
