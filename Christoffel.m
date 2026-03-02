function [vp, vs, pol_vector, Gamma] = Christoffel(n, rho, C)
% Christoffel  Phase velocities and polarization vectors from Christoffel tensor
%
% Inputs
%   n   : propagation direction (3x1 or 1x3), not necessarily normalized
%   rho : density (scalar, >0)
%   C   : stiffness tensor, either
%         - 6x6   (Voigt matrix)
%         - 3x3x3x3 (4th-order tensor)
%
% Outputs
%   vp         : P-wave phase velocity (scalar)
%   vs         : [1x2] S-wave phase velocities (sorted descending, so vs(1)>=vs(2))
%   pol_vector : 3x3 matrix, columns = polarization vectors [P, S1, S2]
%   Gamma      : 3x3 Christoffel matrix
%
% Notes
%   - Eigenvalues are sorted in descending order: lambdaP >= lambdaS1 >= lambdaS2.
%   - For isotropic media, S1 and S2 are degenerate and their polarization basis is not unique.
%   - Assumes the same Voigt convention as your Mat6toTens4() function.

% Input checks
if ~(isscalar(rho) && isnumeric(rho) && isfinite(rho) && rho > 0)
    error('rho must be a positive scalar.');
end

n = double(n(:));
if numel(n) ~= 3
    error('n must be a 3-component vector.');
end
nn = norm(n);
if ~(isfinite(nn) && nn > 0)
    error('n must be nonzero.');
end
n = n / nn;  % normalize direction


% Convert C to 4th-order tensor if needed
if ismatrix(C)
    if ~isequal(size(C), [6 6])
        error('If C is a matrix, it must be 6x6.');
    end
    C4 = MatrixToTensorConversion(C);
elseif ndims(C) == 4
    if ~isequal(size(C), [3 3 3 3])
        error('If C is 4th-order, it must be 3x3x3x3.');
    end
    C4 = C;
else
    error('C must be either 6x6 or 3x3x3x3.');
end

% Christoffel tensor Gamma_ik = C_ijkl n_j n_l
Gamma = zeros(3,3);
for i = 1:3
    for k = 1:3
        s = 0.0;
        for j = 1:3
            for l = 1:3
                s = s + C4(i,j,k,l) * n(j) * n(l);
            end
        end
        Gamma(i,k) = s;
    end
end

% Force symmetry (numerical cleanup)
Gamma = 0.5 * (Gamma + Gamma.');

% Eigen decomposition (exact for 3x3)
[V, D] = eig(Gamma, 'vector');   % D is 3x1 vector of eigenvalues
D = real(D);
V = real(V);

% Sort descending eigenvalues => P first
[Dsort, idx] = sort(D, 'descend');
V = V(:, idx);

% Numerical safeguard against tiny negative values from roundoff
Dsort(Dsort < 0 & Dsort > -1e-12*max(1,max(abs(Dsort)))) = 0;
if any(Dsort < 0)
    warning('Christoffel: negative eigenvalue(s) found. Check C, rho, or conventions.');
end

v = sqrt(max(Dsort,0) / rho);

vp = v(1);
vs = v(2:3).';   % row vector [vs1 vs2]

% Make P polarization roughly aligned with n (sign convention only)
if dot(V(:,1), n) < 0
    V(:,1) = -V(:,1);
end

pol_vector = V;  % columns: [P, S1, S2]

end
