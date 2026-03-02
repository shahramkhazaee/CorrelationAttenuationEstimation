function [sigma,eta] = attenuationCoeffsWeaverDiscretized(Nreal, N, f, Np, Rmax, material, vectors, fnames, homoType, etta)
% Monte Carlo estimate of scattering-induced attenuation coefficients 
% without the assumption of isotropic correlation structure
%
% This function computes the modal differential scattering cross-sections
% (P/S mode conversions) following the formulation used in Norouzian (2019),
% using:
%   1) a Monte Carlo sampling of scattering vectors in a ball of radius Rmax,
%   2) a stationary two-point correlation estimate eta(r,theta,phi),
%   3) a covariance tensor built from grain-wise rotated stiffness tensors.
%
% Inputs
%   Nreal    : vector of realization IDs (e.g., [1 3 7] or 1:10)
%   N        : number of Monte Carlo samples for the scattering integral
%   f        : wave frequency (Hz)
%   Np       : number of point pairs for each pointwise TPCF estimate
%   Rmax     : maximum lag distance used in the correlation/integration (m)
%   material : struct with fields:
%                - elasticConstants : [c11 c12 c44] (cubic, in GPa), or
%                                     [c11 c12 c13 c33 c44] (hexagonal, in GPa)
%                - density          : mass density (kg/m^3)
%   vectors  : struct with fields:
%                - inc_wave_vector  : incident propagation direction (3x1 or 1x3)
%                - sca_wave_vector  : scattered propagation direction (3x1 or 1x3)
%   fnames   : struct with field:
%                - vol_fname : base path/prefix used to read:
%                              [vol_fname num2str(realID) '.stcell']
%                              [vol_fname num2str(realID) '.stvox']
%   etta     : optional precomputed eta array of size [length(Nreal) x N]
%
% Outputs
%   sigma : struct containing the 9 modal cross-sections:
%           sigmaPP, sigmaPS1, sigmaPS2,
%           sigmaS1P, sigmaS1S1, sigmaS1S2,
%           sigmaS2P, sigmaS2S1, sigmaS2S2
%   eta   : estimated (or provided) stationary correlation values used in the MC sum
%
% Notes
%   - If etta is provided, the pointwise TPCF estimation is skipped.
%   - The .stvox file is used to determine the domain bounds and build the
%     voxel-index interpolant for the stationary pointwise TPCF estimator.
%
% Reference
%   Norouzian, M. and Turner, J.A., 2019. Ultrasonic wave propagation 
%   predictions for polycrystalline materials using three-dimensional 
%   synthetic microstructures: Attenuation. JASA, 145(4), pp.2181-2191.
%
% Read inputs
% Material & Parameter Setup
if numel(material.elasticConstants) == 3
    c11 = material.elasticConstants(1);
    c12 = material.elasticConstants(2);
    c44 = material.elasticConstants(3);
    C_local = stiffnessMatrix('cubic', c11, c12, c44) * 1e9;
elseif numel(material.elasticConstants) == 5
    c11 = material.elasticConstants(1);
    c12 = material.elasticConstants(2);
    c13 = material.elasticConstants(3);
    c33 = material.elasticConstants(4);
    c44 = material.elasticConstants(5);
    C_local = stiffnessMatrix('hexagonal', c11, c12, c13, c33, c44) * 1e9;
else
    error('Size of elasticConstants must be either 3 (cubic) or 5 (hexagonal).');
end

rho = material.density;

inc_wave_vector = vectors.inc_wave_vector;
sca_wave_vector = vectors.sca_wave_vector;

vol_fname = fnames.vol_fname;

omega = 2*pi*f;  % angular frequency

if nargin < 9
    homoType = 'voigt'; % default homogenization technique is Voigt
end
homoType = lower(string(homoType));

% initialization of outputs
sigmaPP = zeros(length(Nreal),1); sigmaPS1 = zeros(length(Nreal),1); sigmaPS2 = zeros(length(Nreal),1);
sigmaS1P = zeros(length(Nreal),1); sigmaS1S1 = zeros(length(Nreal),1); sigmaS1S2 = zeros(length(Nreal),1);
sigmaS2P = zeros(length(Nreal),1); sigmaS2S1 = zeros(length(Nreal),1); sigmaS2S2 = zeros(length(Nreal),1);

% check if the TPCF is already available
useProvidedEta = (nargin == 10);
if useProvidedEta
    eta = etta;
    assert(isequal(size(eta), [length(Nreal), N]), 'etta must be size [length(Nreal), N].');
else
    eta = zeros(length(Nreal), N);
end

% loop over realizations
for j=1:length(Nreal)

    rID = Nreal(j);   % actual realization ID

    stcellPath = [vol_fname num2str(rID) '.stcell'];
    fid = fopen(stcellPath, 'r');
    if fid < 0
        error('Could not open file: %s', stcellPath);
    end
    info_Neper = textscan(fid,'%f %*q %f %f %f %*[^\n]');
    fclose(fid);

    vol_list = [info_Neper{1}];

    % number of grains and volume of the sample
    Ng = length(vol_list);
    vol = sum(vol_list); % volumes in Neper are in micron^3 or m^3????

    % Read corresponding Neper voxel file to get true domain extents
    stvoxPath = [vol_fname num2str(rID) '.stvox'];
    vox = read_stvox_domain(stvoxPath, Rmax);

    % Use geometric lengths from .stvox
    bounds = vox.bounds;
    is2D = vox.is2D;
    Rmax_used = vox.Rmax_used;

    % Voxel id list and interpolant (for TPCF estimation)
    cell_id = vox.cell_id;
    interp_function = vox.interp_function;
    
    rng(rID + 1000000, 'twister'); % Reproducibility

    lambda = rand(N,3);

    r_n =  Rmax_used*lambda(:,1).^(1/3);
    mu_n = 2*lambda(:,2)-1;
    phi_n = 2*pi*lambda(:,3);
    s = sqrt(max(0, 1-mu_n.^2));
    ux = s .* cos(phi_n);
    uy = s .* sin(phi_n);
    uz = mu_n;

    if ~useProvidedEta
        % compute eta(j,:)
        % estimate eta for this realization using the voxel-based stationary TPCF estimator
        for i=1:N
            theta_i = acos(mu_n(i)); % 3D polar angle (unused if is2D)
            phi_i   = phi_n(i);      % azimuth
            eta(j,i) = TPCF_stationary_point(Np, r_n(i), theta_i, phi_i, interp_function, cell_id, bounds, is2D);
            if rem(i,500)==0
                disp(['eta:' num2str(i) '/' num2str(N) ' in ' num2str(j) '/' num2str(length(Nreal)) ' realization done!']);
            end
        end
    end

    % Euler angles for all grains of realization j
    % Calculate C0 and C_ave
    rng(rID, 'twister'); % Reproducibility

    phi1 = 2*pi*rand(1,Ng);
    phi = acos(-1+2*rand(1,Ng));
    phi2 = 2*pi*rand(1,Ng);

    R = rotationMatrix(phi1,phi,phi2); % Matrix of all rotations (Ng) for realization j

    C_rotated = stiffnessMatrixLocalToGlobal(C_local, R); % Matrix of all global C (Ng) for realization j

    % Effective medium according to selected homogenization scheme
    C_eff = homogenizedStiffnessDiscrete(C_rotated, vol_list, material, homoType);

    % fluctuation part of the tensor
    deltaC_all = C_rotated - C_eff;

    cov_tensor = zeros(3,3,3,3,3,3,3,3);
    for ig = 1:Ng
        deltaC = MatrixToTensorConversion(deltaC_all(:,:,ig));
        cov_tensor = cov_tensor + vol_list(ig) * tensorOuterProduct(deltaC, deltaC);
    end
    cov_tensor = cov_tensor / vol;

    % calculate incident/scattered polarization vector/phase velocities
    [vp_inc,vs_inc,inc_pol_vector] = Christoffel(inc_wave_vector,rho,C_eff);
    [vp_sca,vs_sca,sca_pol_vector] = Christoffel(sca_wave_vector,rho,C_eff);

    % Eq. (6) of Norouzian 2019, using actual capped radius
    coeff = Rmax_used^3 * omega^4 / (24*pi*rho^2*N);

    % P to P
    coeff_total = coeff/vp_inc^3/vp_sca^5;
    q = (omega/vp_inc)*inc_wave_vector-(omega/vp_sca)*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,1),sca_pol_vector(:,1),inc_wave_vector,sca_wave_vector);
    sigmaPP(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % P to S1
    coeff_total = coeff/vp_inc^3/vs_sca(1)^5;
    q = (omega/vp_inc)*inc_wave_vector-(omega/vs_sca(1))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,1),sca_pol_vector(:,2),inc_wave_vector,sca_wave_vector);
    sigmaPS1(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % P to S2
    coeff_total = coeff/vp_inc^3/vs_sca(2)^5;
    q = (omega/vp_inc)*inc_wave_vector-(omega/vs_sca(2))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,1),sca_pol_vector(:,3),inc_wave_vector,sca_wave_vector);
    sigmaPS2(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S1 to P
    coeff_total = coeff/vs_inc(1)^3/vp_sca^5;
    q = (omega/vs_inc(1))*inc_wave_vector-(omega/vp_sca)*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,2),sca_pol_vector(:,1),inc_wave_vector,sca_wave_vector);
    sigmaS1P(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S1 to S1
    coeff_total = coeff/vs_inc(1)^3/vs_sca(1)^5;
    q = (omega/vs_inc(1))*inc_wave_vector-(omega/vs_sca(1))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,2),sca_pol_vector(:,2),inc_wave_vector,sca_wave_vector);
    sigmaS1S1(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S1 to S2
    coeff_total = coeff/vs_inc(1)^3/vs_sca(2)^5;
    q = (omega/vs_inc(1))*inc_wave_vector-(omega/vs_sca(2))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,2),sca_pol_vector(:,3),inc_wave_vector,sca_wave_vector);
    sigmaS1S2(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S2 to P
    coeff_total = coeff/vs_inc(2)^3/vp_sca^5;
    q = (omega/vs_inc(2))*inc_wave_vector-(omega/vp_sca)*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,3),sca_pol_vector(:,1),inc_wave_vector,sca_wave_vector);
    sigmaS2P(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S2 to S1
    coeff_total = coeff/vs_inc(2)^3/vs_sca(1)^5;
    q = (omega/vs_inc(2))*inc_wave_vector-(omega/vs_sca(1))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,3),sca_pol_vector(:,2),inc_wave_vector,sca_wave_vector);
    sigmaS2S1(j) = coeff_total*lambda_proj*eta(j,:) * phase;
    
    % S2 to S2
    coeff_total = coeff/vs_inc(2)^3/vs_sca(2)^5;
    q = (omega/vs_inc(2))*inc_wave_vector-(omega/vs_sca(2))*sca_wave_vector;
    phase = exp(-1i * r_n .* (ux*q(1) + uy*q(2) + uz*q(3)));
    lambda_proj = projectCovarianceTensor(cov_tensor,inc_pol_vector(:,3),sca_pol_vector(:,3),inc_wave_vector,sca_wave_vector);
    sigmaS2S2(j) = coeff_total*lambda_proj*eta(j,:) * phase;
end

sigma = struct('sigmaPP',sigmaPP,'sigmaPS1',sigmaPS1,'sigmaPS2',sigmaPS2,...
               'sigmaS1P',sigmaS1P,'sigmaS1S1',sigmaS1S1,'sigmaS1S2',sigmaS1S2,...
               'sigmaS2P',sigmaS2P,'sigmaS2S1',sigmaS2S1,'sigmaS2S2',sigmaS2S2);

end

% -------------------------------------------------------------------------
% Helper functions
function vox = read_stvox_domain(stvox_filename, Rmax)
% Read Neper *.stvox (voxel stats) and return domain extents + interpolant.
%
% Outputs (fields):
%   cell_id         [Nv x 1]
%   bounds          [3 x 2] = [xmin xmax; ymin ymax; zmin zmax]
%   is2D            true if z-extent is (numerically) zero
%   Rmax_used       min(Rmax, 0.5*min(domain lengths))
%   interp_function scatteredInterpolant mapping (x,y[,z]) -> voxel index

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

info_vox = readmatrix(stvox_filename, 'FileType', 'text');
if isempty(info_vox) || size(info_vox,2) < 3
    error('File "%s" does not look like a Neper *.stvox with columns: id x y [z].', stvox_filename);
end

cell_id = info_vox(:,1);

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

if is2D
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly]));
    interp_function = scatteredInterpolant(coords(:,1), coords(:,2), ...
        (1:size(coords,1))', 'nearest', 'nearest');
else
    Rmax_used = min(Rmax, 0.5 * min([Lx, Ly, Lz]));
    interp_function = scatteredInterpolant(coords(:,1), coords(:,2), coords(:,3), ...
        (1:size(coords,1))', 'nearest', 'nearest');
end

vox = struct();
vox.cell_id = cell_id;
vox.bounds = [xmin xmax; ymin ymax; zmin zmax];
vox.is2D = is2D;
vox.Rmax_used = Rmax_used;
vox.interp_function = interp_function;
end

% -------------------------------------------------------------------------
function w = TPCF_stationary_point(Np, r, theta, phi, interp_function, cell_id, bounds, is2D)
% TPCF_stationary_point
% Estimate stationary same-grain probability at a single lag/direction.
%
% Inputs
%   Np             : number of random pairs
%   r              : lag distance (same unit as stvox coordinates)
%   theta, phi     : direction angles (theta ignored in 2D)
%   interp_function: scatteredInterpolant returning voxel row indices
%   cell_id        : grain/cell IDs for each voxel row
%   bounds         : [xmin xmax; ymin ymax; zmin zmax]
%   is2D           : logical flag
%
% Output
%   w              : estimated probability that two points belong to the same grain

xmin = bounds(1,1); xmax = bounds(1,2);
ymin = bounds(2,1); ymax = bounds(2,2);
zmin = bounds(3,1); zmax = bounds(3,2);

Lx = xmax - xmin;
Ly = ymax - ymin;
Lz = zmax - zmin;

if is2D
    x1 = [xmin + Lx*rand(Np,1), ymin + Ly*rand(Np,1)];
    dx = [r*cos(phi), r*sin(phi)];
    x2 = x1 + dx;

    in = (x2(:,1) >= xmin) & (x2(:,1) <= xmax) & ...
        (x2(:,2) >= ymin) & (x2(:,2) <= ymax);

    if ~any(in)
        w = NaN;
        return;
    end

    k1 = interp_function(x1(in,1), x1(in,2));
    k2 = interp_function(x2(in,1), x2(in,2));

else
    x1 = [xmin + Lx*rand(Np,1), ymin + Ly*rand(Np,1), zmin + Lz*rand(Np,1)];
    dx = [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)];
    x2 = x1 + dx;

    in = (x2(:,1) >= xmin) & (x2(:,1) <= xmax) & ...
        (x2(:,2) >= ymin) & (x2(:,2) <= ymax) & ...
        (x2(:,3) >= zmin) & (x2(:,3) <= zmax);

    if ~any(in)
        w = NaN;
        return;
    end

    k1 = interp_function(x1(in,1), x1(in,2), x1(in,3));
    k2 = interp_function(x2(in,1), x2(in,2), x2(in,3));
end

% Cleanup indices
ok = ~isnan(k1) & ~isnan(k2);
k1 = k1(ok);
k2 = k2(ok);

if isempty(k1)
    w = NaN;
    return;
end

k1 = round(k1);
k2 = round(k2);

nV = numel(cell_id);
ok = (k1 >= 1) & (k1 <= nV) & (k2 >= 1) & (k2 <= nV);

k1 = k1(ok);
k2 = k2(ok);

if isempty(k1)
    w = NaN;
    return;
end

% compare grain IDs
w = mean(cell_id(k1) == cell_id(k2));
end
