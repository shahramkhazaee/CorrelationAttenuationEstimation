function [w, r, meta] = TPCF_stationary_singleRealization(stvox_filename, Np, Nr, theta, phi, Rmax)
% Directional stationary TPCF estimation for one Neper realization (2D or 3D).
%
%   [w, r] = TPCF_stationary_singleRealization(stvox_filename, Np, Nr)
%   [w, r] = TPCF_stationary_singleRealization(stvox_filename, Np, Nr, theta, phi)
%   [w, r] = TPCF_stationary_singleRealization(stvox_filename, Np, Nr, theta, phi, Rmax)
%   [w, r, meta] = TPCF_stationary_singleRealization(...)
%
% This routine estimates the stationary (translation-invariant) two-point
% correlation along a chosen direction:
%   w(r; u) = E[ 1_{id(x)=id(x + r*u)} ]
% where u is fixed by (theta,phi) in 3D, or by phi in 2D.
%
% Inputs
%   stvox_filename : Neper voxel file name (*.stvox). If no extension is given,
%                   ".stvox" is appended automatically.
%   Np             : number of sampled pairs per lag distance (positive scalar)
%   Nr             : number of lag distances (positive integer scalar)
%   theta, phi     : angles in radians (optional)
%                   - 3D: theta = angle with +z axis, phi = azimuth
%                   - 2D: theta is ignored; phi is the in-plane angle
%                   You can pass scalars (fixed direction) or vectors of length Nr.
%                   If omitted, a random direction is chosen once and reused for all r.
%   Rmax (optional): maximum lag distance. If omitted, Rmax = Inf, then capped by the domain size.
%
% Outputs
%   w    : [Nr x 1] estimated stationary TPCF values
%   r    : [Nr x 1] lag distance vector (linspace from 0 to Rmax_used)
%   meta : struct with fields: is2D, bounds, L, Rmax_used, query (theta,phi)
%
% Notes
%   - 2D vs 3D is detected from the z-range in the *.stvox coordinates.
%   - Uses nearest-neighbor scatteredInterpolant to map (x,y[,z]) -> voxel index.

% Estimate TPCF from Neper outputs for multiple realizations.
%
% Interface:
%   [eta, r, meta] = TPCF(realIDs, Np, Nr, vol_fname)
%   [eta, r, meta] = TPCF(realIDs, Np, Nr, vol_fname, method)
%   [eta, r, meta] = TPCF(realIDs, Np, Nr, vol_fname, method, phase)
%   [eta, r, meta] = TPCF(realIDs, Np, Nr, vol_fname, method, phase, Rmax)
%
% Extra options (mainly for stationary directional cuts):
%   [eta, r, meta] = TPCF(..., 'theta', theta0, 'phi', phi0)
%
% Inputs
%   realIDs   : array of realization indices (e.g. 1:10)
%   Np        : number of random point pairs per evaluation
%   Nr        : number of lag distances
%   vol_fname : base filename prefix (Option A). Example:
%              vol_fname = "2D..._sample" -> reads "..._sample1.stvox", "..._sample2.stvox", ...
%
% Optional inputs
%   method : 'isotropic' (default) or 'stationary'
%   phase  : 'monophase' (default), 'biphase_total', 'biphase_phase1', 'biphase_phase2'
%   Rmax   : maximum lag distance (same unit as Neper coordinates). If omitted, capped by domain size.
%   theta  : (stationary, 3D only) fixed polar angle with +z axis in radians
%   phi    : (stationary, 2D/3D) fixed azimuth angle in radians
%
% Outputs
%   eta  : [nReal x Nr] TPCF values
%   r    : [Nr x 1] lag distance vector (used for both isotropic and stationary)
%   meta : struct array (length nReal) with fields:
%          - files (stvox, stcell when applicable)
%          - bounds, L, is2D
%          - Rmax_cap (0.5*min(L))
%          - Rmax_used (global, possibly capped)
%          - phase, method


% ---- defaults for optional inputs
if nargin < 6
    Rmax = Inf;
end
if nargin < 5
    phi = [];
end
if nargin < 4
    theta = [];
end

% input checks
if ~(isscalar(Np) && isnumeric(Np) && isfinite(Np) && (Np > 0))
    error('Np must be a positive scalar.');
end
if ~(isscalar(Nr) && isnumeric(Nr) && isfinite(Nr) && (Nr > 0) && (Nr == floor(Nr)))
    error('Nr must be a positive integer scalar.');
end

% normalize filename and extension
stvox_filename = strtrim(string(stvox_filename));
[~,~,ext] = fileparts(stvox_filename);

if ext == ""
    stvox_filename = stvox_filename + ".stvox";
elseif ~strcmpi(ext, ".stvox")
    error('The file extension should be .stvox (got: %s).', ext);
end

if ~isfile(stvox_filename)
    error('File not found: %s', stvox_filename);
end

% read voxel file
info_vox = readmatrix(stvox_filename, 'FileType', 'text');
if isempty(info_vox) || size(info_vox,2) < 3
    error('File "%s" does not look like a Neper *.stvox with columns: id x y [z].', stvox_filename);
end

cell_id = info_vox(:,1);

% coordinates: accept both (id,x,y) and (id,x,y,z)
if size(info_vox,2) >= 4
    coords = info_vox(:,2:4);
else
    coords = [info_vox(:,2:3), zeros(size(info_vox,1),1)];
end

xmin = min(coords(:,1)); xmax = max(coords(:,1));
ymin = min(coords(:,2)); ymax = max(coords(:,2));
zmin = min(coords(:,3)); zmax = max(coords(:,3));

Lx = xmax - xmin;
Ly = ymax - ymin;
Lz = zmax - zmin;

tolZ = 1e-12 * max([Lx, Ly, 1]);
is2D = (Lz <= tolZ);

% cap Rmax to limit boundary rejections
if is2D
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly]));
    % Build interpolant: (x,y) -> voxel index
    ns = createns(coords(:,1:2), 'NSMethod','kdtree');
else
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly, Lz]));
    % Build interpolant: (x,y,z) -> voxel index
    ns = createns(coords, 'NSMethod','kdtree');
end

% nearest voxel index
interp_function = @(q) knnsearch(ns, q);

% lag distances
r = linspace(0, Rmax_used, Nr).';
w = NaN(Nr,1);

% build query angles (fixed direction unless user provides vectors)
if isempty(phi)
    % choose one random direction for the whole curve
    if is2D
        phi = 2*pi*rand * ones(Nr,1);
        theta = NaN(Nr,1);
    else
        phi = 2*pi*rand * ones(Nr,1);
        theta = acos(2*rand - 1) * ones(Nr,1); % cos(theta) uniform
    end
else
    % phi provided (scalar or vector)
    if isscalar(phi)
        phi = double(phi) * ones(Nr,1);
    else
        phi = double(phi(:));
        if numel(phi) ~= Nr
            error('phi must be a scalar or a vector of length Nr.');
        end
    end

    if is2D
        theta = NaN(Nr,1);
    else
        if isempty(theta)
            error('For 3D, theta must be provided (scalar or vector of length Nr).');
        end
        if isscalar(theta)
            theta = double(theta) * ones(Nr,1);
        else
            theta = double(theta(:));
            if numel(theta) ~= Nr
                error('theta must be a scalar or a vector of length Nr.');
            end
        end
    end
end

% main loop over lag distances
for i = 1:Nr
    rr = r(i);

    if is2D
        % x1 uniform in rectangle, z fixed
        x1 = [xmin + Lx*rand(Np,1), ymin + Ly*rand(Np,1)];
        dx = [rr*cos(phi(i)), rr*sin(phi(i))];
        x2 = x1 + dx;

        in = (x2(:,1) >= xmin) & (x2(:,1) <= xmax) & ...
             (x2(:,2) >= ymin) & (x2(:,2) <= ymax);

        x1v = x1(in,:);
        x2v = x2(in,:);

        if isempty(x1v)
            w(i) = NaN;
            continue;
        end

        k1 = interp_function(x1v);
        k2 = interp_function(x2v);

    else
        % x1 uniform in box
        x1 = [xmin + Lx*rand(Np,1), ymin + Ly*rand(Np,1), zmin + Lz*rand(Np,1)];
        dx = [rr*sin(theta(i))*cos(phi(i)), rr*sin(theta(i))*sin(phi(i)), rr*cos(theta(i))];
        x2 = x1 + dx;

        in = (x2(:,1) >= xmin) & (x2(:,1) <= xmax) & ...
             (x2(:,2) >= ymin) & (x2(:,2) <= ymax) & ...
             (x2(:,3) >= zmin) & (x2(:,3) <= zmax);

        x1v = x1(in,:);
        x2v = x2(in,:);

        if isempty(x1v)
            w(i) = NaN;
            continue;
        end

        k1 = interp_function(x1v);
        k2 = interp_function(x2v);
    end

    % Remove NaNs (should be rare with nearest extrapolation)
    ok = ~isnan(k1) & ~isnan(k2);
    k1 = k1(ok); k2 = k2(ok);
    if isempty(k1)
        w(i) = NaN;
        continue;
    end

    % Ensure indices are valid
    k1 = round(k1);
    k2 = round(k2);
    nV = numel(cell_id);
    ok = (k1 >= 1) & (k1 <= nV) & (k2 >= 1) & (k2 <= nV);
    k1 = k1(ok); k2 = k2(ok);

    if isempty(k1)
        w(i) = NaN;
        continue;
    end

    w(i) = mean(cell_id(k2) == cell_id(k1));

    if rem(i,10) == 0
        disp(['calc:' num2str(i) '/' num2str(Nr) ' lag distances done!']);
    end
end

if nargout >= 3
    meta = struct();
    meta.is2D = is2D;
    meta.bounds = [xmin xmax; ymin ymax; zmin zmax];
    meta.L = [Lx Ly Lz];
    meta.Rmax_used = Rmax_used;
    meta.query = struct('theta', theta, 'phi', phi);
end

end
