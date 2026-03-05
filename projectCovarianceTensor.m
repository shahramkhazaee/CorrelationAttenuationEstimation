function out = projectCovarianceTensor(cov_tensor, inc_pol_vector, sca_pol_vector, inc_wave_vector, sca_wave_vector)
% Projects an 8th-order covariance tensor onto incident/scattered polarization
% and wave vectors.
%
% Required functions
%   - tensorOuterProduct
%
% This computes the scalar:
%   out = C_{ijklmnop} * u_i * p_j * s_k * v_l * u_m * p_n * s_o * v_p
% where:
%   u = incident polarization vector
%   p = incident wave vector (direction)
%   s = scattered wave vector (direction)
%   v = scattered polarization vector
%
% Inputs
%   cov_tensor      : 3x3x3x3x3x3x3x3 covariance tensor
%   inc_pol_vector  : incident polarization vector (3x1 or 1x3)
%   sca_pol_vector  : scattered polarization vector (3x1 or 1x3)
%   inc_wave_vector : incident wave vector (3x1 or 1x3), usually normalized
%   sca_wave_vector : scattered wave vector (3x1 or 1x3), usually normalized
%
% Output
%   out : scalar projection value

% Input checks
if ndims(cov_tensor) ~= 8 || any(size(cov_tensor) ~= 3)
    error('cov_tensor must be a 3x3x3x3x3x3x3x3 array.');
end

u = inc_pol_vector(:);
v = sca_pol_vector(:);
p = inc_wave_vector(:);
s = sca_wave_vector(:);

if numel(u) ~= 3 || numel(v) ~= 3 || numel(p) ~= 3 || numel(s) ~= 3
    error('All input vectors must have exactly 3 components.');
end

% Build X(i,j,k,l) = u_i * p_j * s_k * v_l
X = reshape(u, [3 1 1 1]) .* ...
    reshape(p, [1 3 1 1]) .* ...
    reshape(s, [1 1 3 1]) .* ...
    reshape(v, [1 1 1 3]);

% Tensor outer product
XX = tensorOuterProduct(X, X);   % size 3x3x3x3x3x3x3x3

% Contractions
out = sum(cov_tensor(:) .* XX(:));

end
