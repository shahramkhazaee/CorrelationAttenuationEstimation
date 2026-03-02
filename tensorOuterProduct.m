function C = tensorOuterProduct(A, B)
% Tensor outer product (without contraction)
%
% Returns an array C of size [size(A), size(B)] such that
%   C(i1,...,im,j1,...,jn) = A(i1,...,im) * B(j1,...,jn)
%
% Notes
%   - Size of the output is [size(A), size(B)].
%   - Compatible with MATLAB R2016b+

    Aresh = reshape(A, [size(A), ones(1, ndims(B))]);
    Bresh = reshape(B, [ones(1, ndims(A)), size(B)]);

    C = Aresh .* Bresh;   % implicit expansion
end
