function out = attenuationCoefficientsWeaverMonodisperse(material, source)
% Ultrasonic attenuation calculation using Weaver's model (Born Approx) for
% 3D polycrystals without taking into account the grain size distribution
%
% Note: Requires external functions: 
%   - stiffnessMatrix.m
%   - phaseVelocityHomogenized.m
%   - InnerProducts.m
%
%  INPUTS
%  material (struct variable)
%      elasticConstants [c11 c12 c44 ...] : elastic constant (GPa)
%      density                            : material density (kg/m^3)
%      meanD                              : vector of (mean) grain sizes
%      PSDF                               : power spectral density function
%                                          'exp','gaussian','spherical','sk'
%  source (struct variable)
%      Pwavenumber                        : vector of wavenumber of P waves
%
%  OUTPUTS
%  attenuation coefficients in Np/m: 
%  alphaPP  P-to-P wave attenuation coefficient
%  alphaPS  P-to-S wave attenuation coefficient
%  alphaSP  S-to-P wave attenuation coefficient
%  alphaSS  S-to-S wave attenuation coefficient
%
%  Note: 
%  alphaP = alphaPP + alphaPS;
%  alphaS = alphaSS + alphaSP;

% References: 
% 1. Kube & Turner 2015, Acoustic attenuation coefficients for polycrystalline 
% materials containing crystallites of any symmetry class
% 2. Weaver 1990, Diffusivity of ultrasound in polycrystals
% 3. Khazaie et al, 2016, Influence of the spatial correlation structure of 
% an elastic random medium on its scattering properties

% Example:
% c11 = 103.4; c12 = 57.1; c44 = 28.6; rho = 2700; % Al Stanke Kino
% material.elasticConstants = [c11 c12 c44]; material.density = rho;
% material.meanD = linspace(10e-6,100e-6,20);
% material.PSDF = 'spherical';
% freqs = linspace(1e6,5e6,30);
% V = phaseVelocityHomogenized(stiffnessMatrix('cubic',c11,c12,c44)*1e9,rho,'Voigt');
% source.Pwavenumber = 2*pi*freqs/V(1);
% out = attenuationCoefficientsWeaverMonodisperse(material, source);
% figure; plot(source.Pwavenumber, out.alphaPP + out.alphaPS)

% Material & Parameter Setup
if numel(material.elasticConstants) == 3
    c11 = material.elasticConstants(1);
    c12 = material.elasticConstants(2);
    c44 = material.elasticConstants(3);
    Cloc = stiffnessMatrix('cubic', c11, c12, c44) * 1e9;
elseif numel(material.elasticConstants) == 5
    c11 = material.elasticConstants(1);
    c12 = material.elasticConstants(2);
    c13 = material.elasticConstants(3);
    c33 = material.elasticConstants(4);
    c44 = material.elasticConstants(5);
    Cloc = stiffnessMatrix('hexagonal', c11, c12, c13, c33, c44) * 1e9;
else
    error('elasticConstants must be length 3 (cubic) or 5 (hexagonal).');
end

rho = material.density;
diameters = material.meanD(:); % Column vector of grain sizes
TPCF_type = material.PSDF;
kp = source.Pwavenumber(:).'; % Row vector of wavenumbers

% Initialize Outputs
% Size: [Num_GrainSizes x Num_Wavenumbers]
Nd = length(diameters);
Nk = length(kp);
out.alphaPP = zeros(Nd, Nk);
out.alphaPS = zeros(Nd, Nk);
out.alphaSP = zeros(Nd, Nk);
out.alphaSS = zeros(Nd, Nk);

% Velocities & Inner Products
V0 = phaseVelocityHomogenized(Cloc, rho, 'Voigt');
Vp0 = V0(1);
Vs0 = V0(2);
K = Vp0 / Vs0;
ks = kp * K; 

[ll, mm, nn] = InnerProducts(Cloc);

% Scattering geometric factors (Kube & Turner 2015, Eq 4)
L = @(x) ll.L0 + ll.L1 * x.^2 + ll.L2 * x.^4;
M = @(x) mm.M0 + mm.M1 * x.^2;
N = @(x) nn.N0 + nn.N1 * x.^2;

% TPCF/PSDF Definitions
% PSDF(z) is the dimensionless Spectral Density where z = k*D_mean
if strcmpi(TPCF_type,'exp')

    PSDF = @(z) 1./(8*pi^2) ./ (1 + z.^2/4).^2;

elseif strcmpi(TPCF_type,'gaussian')

    PSDF = @(z) 1/(8*pi^3) * exp(-z.^2 / (4*pi));

elseif strcmpi(TPCF_type,'spherical') || strcmpi(TPCF_type,'sha') ...
                                      || strcmpi(TPCF_type,'sk')

    PSDF = @(z) psdfSpherical(z);

else
    error(['PSDF type "' TPCF_type '" not implemented in monodisperse routine.']);
end

% Integration loop
abstol = 1e-10; reltol = 1e-6;

% Note: We loop over grain sizes 'd'
for i = 1:Nd
    d = diameters(i);
       
    % PP Scattering (LL)
    % Argument: kp*d * sqrt(2*(1-x))
    funLL = @(x) PSDF(kp .* d .* sqrt(2*(1-x))) .* L(x);
    int_LL = integral(funLL, -1, 1, 'AbsTol', abstol, 'RelTol', reltol, 'ArrayValued', true);
    
    % PS Scattering (LT)
    % Argument: sqrt( (kd)^2 + (ksd)^2 - 2(kd)(ksd)x )
    % Note: ks = kp*K
    kpd = kp * d;
    ksd = ks * d;
    term1 = kpd.^2 + ksd.^2;
    term2 = 2 .* kpd .* ksd;
    
    funLT = @(x) PSDF(sqrt(max(term1 - term2.*x,0))) .* (M(x) - L(x));
    int_LT = integral(funLT, -1, 1, 'AbsTol', abstol, 'RelTol', reltol, 'ArrayValued', true);
    
    % SS Scattering (TT)
    funTT = @(x) PSDF(ksd .* sqrt(2*(1-x))) .* (N(x) - 2*M(x) + L(x));
    int_TT = integral(funTT, -1, 1, 'AbsTol', abstol, 'RelTol', reltol, 'ArrayValued', true);
    
    % Factors
    fPP = (d^3 * pi^2 * kp.^4) / (2 * Vp0^4 * rho^2);
    fPS = (d^3 * pi^2 * ks.^4) / (2 * Vp0^3 * Vs0 * rho^2);
    fSS = (d^3 * pi^2 * ks.^4) / (4 * Vs0^4 * rho^2);
    
    % Weaver's attenuation coefficients (Eqs. 3.35)
    out.alphaPP(i, :) = fPP .* int_LL;
    out.alphaPS(i, :) = fPS .* int_LT;
    out.alphaSP(i, :) = (0.5 / K^2) * out.alphaPS(i, :);
    out.alphaSS(i, :) = fSS .* int_TT;
end

end

% Helper function for the PSDF of the spherical case
function y = psdfSpherical(z)
    y = zeros(size(z));
    small = abs(z) < 1e-8;

    zs = z(~small);
    y(~small) = 3 * (-zs.*cos(zs/2) + 2*sin(zs/2)).^2 ./ (pi^2 * zs.^6);

    y(small) = 1/(48*pi^2);
end
