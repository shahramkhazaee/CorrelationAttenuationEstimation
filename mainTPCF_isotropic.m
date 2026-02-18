function [eta, r, meta] = mainTPCF_isotropic(realIDs, Np, Nr, varargin)
% Estimates isotropic two-point correlation function (TPCF) for multiple
% Neper realizations (2D or 3D), by looping over *.stvox files.
%
% Supported call patterns
%   [eta,r]       = mainTPCF_isotropic(realIDs, Np, Nr, vol_fname)
%   [eta,r]       = mainTPCF_isotropic(realIDs, Np, Nr, vol_fname, Rmax)
%   [eta,r,meta]  = mainTPCF_isotropic(...)
%
% Inputs
%   realIDs   : vector of realization indices (e.g. 1:10)
%   Np        : number of random point pairs per lag distance
%   Nr        : number of lag distances
%   vol_fname : stvox file name prefix (without realization number and extension)
%               Example: "..._sample" for files "..._sample1.stvox", "..._sample2.stvox", ...
%   Rmax      : (optional) requested max lag distance (same unit as Neper coordinates)
%               If omitted, Rmax is computed automatically from the first realization.
% Example
% [eta, r, meta] = mainTPCF_isotropic(1:3, 5e3, 100, '2Dmonomodal_4000x1500_Dbar75_deltaD20_sample',1000);

% ---------------------------
% Parse inputs (Rmax optional)
% ---------------------------
if isempty(varargin)
    error('Missing input "vol_fname".');
end

% Only vol_fname and optional Rmax are allowed
if numel(varargin) > 2
    error('Too many optional inputs. Expected: vol_fname, optionally Rmax.');
end

% vol_fname
if ischar(varargin{1}) || isstring(varargin{1})
    vol_fname = strtrim(string(varargin{1}));
else
    error('Invalid argument. Expected vol_fname as a char or string.');
end

% Optional: Rmax
hasRmax = false;
if numel(varargin) == 2
    if isempty(varargin{2})
        % user explicitly passed [] as Rmax, treat as "not provided"
        hasRmax = false;
    elseif isnumeric(varargin{2}) && isscalar(varargin{2})
        Rmax = varargin{2};
        hasRmax = true;
    else
        error('Invalid optional argument. Expected Rmax as a numeric scalar.');
    end
end

realIDs   = realIDs(:).'; % row
nReal     = numel(realIDs);

if nReal == 0
    error('realIDs is empty. Provide at least one realization index.');
end

eta  = zeros(nReal, Nr);

% Helper to build the base filename (without extension)
make_base = @(rid) vol_fname + string(rid);

% ---------------------------
% First realization
% ---------------------------
rid0 = realIDs(1); % id of the first realization
if hasRmax
    [eta0, r, meta0] = TPCF_isotropic_singleRealization(make_base(rid0), Np, Nr, Rmax);
else
    [eta0, r, meta0] = TPCF_isotropic_singleRealization(make_base(rid0), Np, Nr);
end

eta(1,:) = eta0(:).';
meta = repmat(meta0, 1, nReal);

% Use the first realization to define a global requested Rmax for the rest
Rmax_global = meta0.Rmax_used;
is2D_ref    = meta0.is2D;

% ---------------------------
% Remaining realizations
% ---------------------------
for j = 2:nReal
    rid = realIDs(j);

    [etaj, ~, metaj] = TPCF_isotropic_singleRealization(make_base(rid), Np, Nr, Rmax_global);

    % Do not mix 2D and 3D runs!
    if metaj.is2D ~= is2D_ref
        error('Mixed dimensionality detected across realizations (2D and 3D). Check your input files.');
    end

    eta(j,:) = etaj(:).';
    meta(j)  = metaj;

    if rem(j, 2) == 0 || j == nReal
        disp(['calc:' num2str(j) '/' num2str(nReal) ' realizations done!']);
    end
end

end
