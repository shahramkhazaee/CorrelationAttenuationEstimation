function out = attenuationCoefficientsWeaverPolydisperse(material, source)
% Ultrasonic attenuation calculation using Weaver's model (Born Approx) for
% 3D polycrystals without taking into account the grain size distribution
% Wraps attenuationCoefficientsWeaverMonodisperse.m by integrating over a
% PDF of grain sizes
%
% Supported PDFs: 'lognormal', 'normal', 'gamma', 'weibull'
%
%  INPUTS
%  material (struct variable)
%      elasticConstants [c11 c12 c44 ...] : elastic constant in GPa
%      density                            : material density (kg/m^3)
%      meanD                              : mean of the grain sizes D
%      stdD                               : standard deviation of D
%      PDF                                : PDF of D
%      TPCF                               : TPCF: 'exp','spherical','analytical'
%  source (struct variable)
%      kp                                 : wavenumber of P waves
%
% OUTPUTS
%  attenuation coefficients in Np/m:
%  alphaPP  P-to-P wave attenuation coefficient
%  alphaPS  P-to-S wave attenuation coefficient
%  alphaSP  S-to-P wave attenuation coefficient
%  alphaSS  S-to-S wave attenuation coefficient
%  Note:
%  alphaP = alphaPP + alphaPS;
%  alphaS = alphaSS + alphaSP;

% EXAMPLE
% c11 = 103.4e9; c12 = 57.1e9; c44 = 28.6e9; rho = 2700; % Alu Stanke Kino
% m = 50e-6; delta = 0.5; s = m*delta;
% kp = logspace(log10(1e-2/m), log10(20/m), 50);
% source = struct('Pwavenumber', kp);
% material = struct('elasticConstants', [c11, c12, c44], ...
%                 'density', rho, ...
%                 'PSDF',   'spherical', ...
%                 'meanD',   m, ...
%                 'stdD',    s, ...
%                 'PDF',    'lognormal');
% out = attenuationCoefficientsWeaverPolydisperse(material,source);
% figure; hold on; 
% loglog(kp*m, out.alphaPP + out.alphaPS);
% set(gca,'XScale','log','YScale','log');

% Extract parameters
mD = material.meanD;
if isfield(material, 'stdD')
    sD = material.stdD;
else
    sD = 0;
end

TPCF_type = lower(strtrim(char(material.PSDF)));
PDF_type  = material.PDF;

% Particular case of single grain size (delta) distribution
if sD == 0
    out = attenuationCoefficientsWeaverMonodisperse(material, source);
    return
end

% Integration over distribution of grain sizes
% define a range that covers 99.9% of the PDF
prob_coverage = [0.0005, 0.9995];
[PDF, bounds] = makeGrainDiameterPDF(PDF_type, mD, sD, prob_coverage);

% Discretize grain size Domain
Na = 2e3; % this could be modified to check the convergence
a = linspace(bounds(1), bounds(2), Na);
a = a(:); % ensure column vector

% Estimate attenuation coefficients for all grain sizes in vector 'a'
material.meanD = a;
attcoeffs = attenuationCoefficientsWeaverMonodisperse(material, source);

% Integrate : Integral( alpha(a) * PDF(a) ) da
% We use trapezoidal integration and normalize by the area of PDF on this
% range to be strictly correct
switch TPCF_type
    case 'sk'
        % E[D^3]
        ED3 = integral(@(D) D.^3 .* PDF(D), bounds(1), bounds(2));
        % Volume weighted PDF to be used in the SK model
        pdf_vals = a.^3 .* PDF(a) / ED3;
    otherwise
        pdf_vals = PDF(a);
end
norm_factor = trapz(a, pdf_vals); % Should be close to 1

if round(abs(norm_factor-1),2)>0.02
    fprintf(['Warning: PDF should be refined more! Please increase the value ' ...
        'of Na in the code! \n']);
end

fields = fieldnames(attcoeffs);
out = struct();

for j = 1:length(fields)
    fieldName = fields{j};
    alpha_a_matrix = attcoeffs.(fieldName); % Size: [length(a) x length(kp)]

    % We integrate along dimension 1 (grain sizes)
    % Result is 1 x length(kp)
    out.(fieldName) = trapz(a, alpha_a_matrix .* pdf_vals, 1) / norm_factor;
end

end