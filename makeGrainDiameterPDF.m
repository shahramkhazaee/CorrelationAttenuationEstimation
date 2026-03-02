function varargout = makeGrainDiameterPDF(PDF_type, mD, sD, prob_coverage)
% returns grain-diameter PDF handle and optional bounds
%
% Usage:
%   PDF = makeGrainDiameterPDF(PDF_type, mD, sD)
%   PDF = makeGrainDiameterPDF(PDF_type, mD, sD, prob_coverage)
%   [PDF, bounds] = makeGrainDiameterPDF(PDF_type, mD, sD)
%   [PDF, bounds] = makeGrainDiameterPDF(PDF_type, mD, sD, prob_coverage)
%
% Inputs
%   PDF_type      : 'lognormal', 'truncated_lognormal', 'normal', 
%                   'truncated_normal', 'gamma', 'weibull'
%   mD            : mean diameter (must be > 0)
%   sD            : std diameter  (must be > 0); the zero case corresponds
%   to a Dirac delta distribution around the mean value
%   prob_coverage : optional 1x2 vector like [0.0005 0.9995]
%
% Outputs
%   PDF           : function handle p_D(t)
%   bounds        : optional 1x2 bounds covering prob_coverage
%
% Notes
% For truncated distributions, the bounds are hard-coded.
% The user should change them based on the case study.

% Defaults
if nargin < 4 || isempty(prob_coverage)
    prob_coverage = [0.0005, 0.9995];
end

% Basic input checks
if ~(isscalar(mD) && isfinite(mD) && mD > 0)
    error('mD must be a positive finite scalar.');
end
if ~(isscalar(sD) && isfinite(sD) && sD >= 0)
    error('sD must be a nonnegative finite scalar.');
end
if ~(isnumeric(prob_coverage) && numel(prob_coverage) == 2 && all(isfinite(prob_coverage)))
    error('prob_coverage must be a 1x2 finite numeric vector.');
end
prob_coverage = prob_coverage(:).';
if ~(prob_coverage(1) > 0 && prob_coverage(2) < 1 && prob_coverage(2) > prob_coverage(1))
    error('prob_coverage must satisfy 0 < p1 < p2 < 1.');
end

PDF_type = lower(string(PDF_type));

% coefficient of variation of grain diameters
tol_cv = 1e-6; % tolerance
deltaD = sD / mD;
if abs(deltaD) < tol_cv
    % PDF is a delta function around the mean value
    error('Delta-like PDF (sD/mD < %g): no need to calculate the PDF.', tol_cv);
else
    switch PDF_type
        case 'lognormal'

            sigmalnD = sqrt(log(1 + deltaD^2));
            mulnD = log(mD) - 0.5 * sigmalnD^2;

            bounds = logninv(prob_coverage, mulnD, sigmalnD);
            PDF = @(x) lognpdf(x, mulnD, sigmalnD);

        case {'normal', 'gaussian'}

            bounds = norminv(prob_coverage, mD, sD);
            bounds(bounds < eps) = eps; % physical truncation
            PDF = @(x) normpdf(x, mD, sD);

        case 'gamma'

            theta = (sD^2) / mD;
            k = mD / theta;
            bounds = gaminv(prob_coverage, k, theta);
            PDF = @(x) gampdf(x, k, theta);

        case 'weibull'

            cv_fun = @(k) (sqrt(gamma(1+2./k) - gamma(1+1./k).^2) ./ gamma(1+1./k)) - deltaD;
            k_guess = 1.27 / max(deltaD, eps);

            try
                k = fzero(cv_fun, k_guess);
            catch
                k = k_guess;
            end
            
            lambda = mD / gamma(1 + 1/k);
            bounds = wblinv(prob_coverage, lambda, k);
            PDF = @(x) wblpdf(x, lambda, k);

        case {'truncated_normal','truncated_gaussian','truncatednormal', ...
                                                    'truncatedgaussian'}

            % lower and upper bounds for the truncated distribution 
            % to be modified based on the case study
            a = eps;
            b = norminv(prob_coverage(2), mD, sD); % instead of Inf

            % if b becomes Inf or NaN, set it to 5 standard deviations
            if isinf(b) || isnan(b)
                b = mD + 5*sD; 
            end
            
            % normalization constant
            Z = normcdf(b,mD,sD) - normcdf(a,mD,sD);
            
            if Z <= 0
                error('Invalid truncation interval for truncated normal.');
            end
            
            bounds = [a b];
            PDF = @(x) (normpdf(x,mD,sD)./Z) .* (x>=a & x<=b);

        case {'truncated_lognormal','truncatedlognormal'}

            sigmalnD = sqrt(log(1 + deltaD^2));    
            mulnD = log(mD) - 0.5 * sigmalnD^2;
            
            % lower and upper bounds for the truncated distribution 
            % to be modified based on the case study
            a = eps;
            b = logninv(prob_coverage(2), mD, sD); % instead of Inf

            % if b becomes Inf or NaN, set it to 5 standard deviations
            if isinf(b) || isnan(b)
                b = mulnD + 5*sigmalnD; 
            end
            
            Z = logncdf(b,mulnD,sigmalnD) - logncdf(a,mulnD,sigmalnD);
            
            if Z <= 0
                error('Invalid truncation interval for truncated lognormal.');
            end
            
            bounds = [a b];
            PDF = @(x) (lognpdf(x,mulnD,sigmalnD)./Z) .* (x>=a & x<=b);

        otherwise

            error('PDF type "%s" is not supported. Use lognormal, normal, gamma, or weibull.', PDF_type);
    end
end

% Outputs: 1 output -> only PDF, 2 outputs -> [PDF, bounds]
if nargout <= 1
    varargout = {PDF};
else
    varargout = {PDF, bounds};
end

end
