function [eta, r, meta] = TPCF_isotropic_singleRealization(stvox_filename, Np, Nr, Rmax)
% Estimates isotropic TPCF for one single Neper realization (2D or 3D).
%
%   [eta, r] = TPCF_isotropic_singleRealization(stvox_file, Np, Nr, Rmax)
%   [eta, r, meta] = TPCF_isotropic_singleRealization(...)
%
% Inputs
%   stvox_filename: filename of Neper voxel file (*.stvox) for a single realization
%   Np            : number of random point pairs per lag distance
%   Nr            : number of lag distances
%   Rmax(optional): maximum lag distance (same unit as Neper coordinates)
%
% Outputs
%   eta  : [1 x Nr] estimated TPCF values
%   r    : [1 x Nr] lag distance vector (same unit as Neper coordinates)
%   meta : struct with fields:
%          - is2D
%          - bounds = [xmin xmax; ymin ymax] or [xmin xmax; ymin ymax; zmin zmax]
%          - L      = [Lx Ly] or [Lx Ly Lz]
%          - Rmax_used
%
% Notes
%   - Automatically detects 2D vs 3D from the ranges in the *.stvox file.
%   - Uses scatteredInterpolant with nearest-neighbor interpolation and extrapolation.
%
% Example
% [eta, r, meta] = TPCF_isotropic_singleRealization('100_35_lognormal.stvox', 1e4, 100);

if nargin < 4
    Rmax = Inf;
end

% Check extension of the input file name (from Neper)
stvox_filename = strtrim(string(stvox_filename));
[~,~,ext] = fileparts(stvox_filename);

if ext == ""
    stvox_filename = stvox_filename + ".stvox";
elseif ~strcmpi(ext, ".stvox")
    error("The file extension should be .stvox (got: %s).", ext);
end

if ~isfile(stvox_filename)
    error("File not found: %s", stvox_filename);
end

% Read Neper voxel file
info_vox = readmatrix(stvox_filename, 'FileType', 'text');

if isempty(info_vox) || size(info_vox,2) < 3
    error('File "%s" does not look like a Neper *.stvox with columns: id x y [z].', stvox_filename);
end

cell_id = info_vox(:, 1);

% Coordinates: accept both 2D (id,x,y) and 3D (id,x,y,z). 
% If more columns exist, keep 2:4.
if size(info_vox,2) >= 4
    coords = info_vox(:, 2:4);
end

xmin = min(coords(:,1)); xmax = max(coords(:,1));
ymin = min(coords(:,2)); ymax = max(coords(:,2));
zmin = min(coords(:,3)); zmax = max(coords(:,3));

Lx = xmax - xmin;
Ly = ymax - ymin;
Lz = zmax - zmin;

% Robust 2D detection (Neper 2D has all z = 0)
tolZ = 1e-12 * max([Lx, Ly, 1]);
is2D = (Lz <= tolZ);

% Clamp Rmax to reduce boundary rejections (same spirit as your 2D main function)
if is2D
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly]));
    % Build 2D interpolant: (x,y) -> voxel index
    interp_function = scatteredInterpolant(coords(:,1), coords(:,2), ...
                                (1:size(coords,1))', 'nearest', 'nearest');
else
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly, Lz]));
    % Build 3D interpolant: (x,y,z) -> voxel index
    interp_function = scatteredInterpolant(coords(:,1), coords(:,2), coords(:,3), ...
                                (1:size(coords,1))', 'nearest', 'nearest');
end

% Vector of lag distances and estimated TPCF values
r = linspace(0, Rmax_used, Nr);
eta = zeros(1, Nr);

for i = 1:Nr
    rr = r(i);

    if is2D
        % Reference points uniformly in the rectangle
        x1 = [xmin + rand(Np,1)*Lx, ymin + rand(Np,1)*Ly];

        ang = 2*pi*rand(Np,1);
        x2 = x1 + rr * [cos(ang), sin(ang)];

        % check if inside the domain
        ind = x2(:,1) >= xmin & x2(:,1) <= xmax & ...
              x2(:,2) >= ymin & x2(:,2) <= ymax;

        x1v = x1(ind,:);
        x2v = x2(ind,:);

        if isempty(x1v)
            eta(i) = NaN;
            continue;
        end

        k1 = interp_function(x1v(:,1), x1v(:,2));
        k2 = interp_function(x2v(:,1), x2v(:,2));

    else
        % Reference points uniformly in the box
        x1 = [xmin + rand(Np,1)*Lx, ymin + rand(Np,1)*Ly, zmin + rand(Np,1)*Lz];

        % Uniform random direction on unit sphere
        phi   = 2*pi*rand(Np,1);
        theta = acos(-1 + 2*rand(Np,1));

        dispVec = rr * [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
        x2 = x1 + dispVec;

        ind = x2(:,1) >= xmin & x2(:,1) <= xmax & ...
              x2(:,2) >= ymin & x2(:,2) <= ymax & ...
              x2(:,3) >= zmin & x2(:,3) <= zmax;

        x1v = x1(ind,:);
        x2v = x2(ind,:);

        if isempty(x1v)
            eta(i) = NaN;
            continue;
        end

        k1 = interp_function(x1v(:,1), x1v(:,2), x1v(:,3));
        k2 = interp_function(x2v(:,1), x2v(:,2), x2v(:,3));
    end

    eta(i) = mean(cell_id(k2) == cell_id(k1));

    if rem(i, 10) == 0
        disp(['calc:' num2str(i) '/' num2str(Nr) ' lag distances done!']);
    end
end

if nargout >= 3
    meta = struct();
    meta.is2D = is2D;
    meta.bounds = [xmin xmax; ymin ymax; zmin zmax];
    meta.L = [Lx Ly Lz];
    meta.Rmax_used = Rmax_used;
end

end
