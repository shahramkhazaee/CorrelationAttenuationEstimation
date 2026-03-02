function C = TensorToMatrixConversion(T, type)
% TensorToMatrixConversion Convert a 3x3x3x3 elasticity tensor to a 6x6 matrix
%
%   C = TensorToMatrixConversion(T)
%   C = TensorToMatrixConversion(T, type)
%
% Inputs
%   T    : 3x3x3x3 tensor
%   type : optional string
%          - omitted / '' : direct mapping with index order [11 22 33 23 13 12]
%          - 'Kelvin'      : applies the legacy normalized-shear scaling used
%                           in the companion MatrixToTensorConversion function
%                           (sqrt(2) and 2 factors)
%
% Output
%   C : 6x6 matrix
%
% Notes
%   - Assumes index order [11, 22, 33, 23, 13, 12].
%   - If T does not satisfy minor symmetries, this function picks the
%     components corresponding to (11,22,33,23,13,12) exactly as indexed.

% Input checks
if ndims(T) ~= 4 || any(size(T) ~= [3 3 3 3])
    error('T must be a 3x3x3x3 tensor.');
end

if nargin < 2
    mode = 'direct';
else
    if ~(ischar(type) || isstring(type))
        error('type must be a character vector or string.');
    end
    type = char(string(type));
    if strcmpi(type, 'Kelvin')
        mode = 'scaled';   % backward-compatible behavior
    else
        error('Unknown type "%s". Use ''Voigt'' or omit the argument.', type);
    end
end

% Pair list in Voigt order [11 22 33 23 13 12]
pairs = [1 1;
         2 2;
         3 3;
         2 3;
         1 3;
         1 2];

C = zeros(6,6);

for I = 1:6
    i = pairs(I,1);
    j = pairs(I,2);

    for J = 1:6
        k = pairs(J,1);
        l = pairs(J,2);

        val = T(i,j,k,l);

        if strcmp(mode, 'direct')
            C(I,J) = val;
        else
            % Inverse of MatrixToTensorConversion(...,'Voigt') scaling
            if I <= 3 && J <= 3
                C(I,J) = val;
            elseif (I <= 3 && J >= 4) || (I >= 4 && J <= 3)
                C(I,J) = sqrt(2) * val;
            else % I>=4 && J>=4
                C(I,J) = 2 * val;
            end
        end
    end
end

end
