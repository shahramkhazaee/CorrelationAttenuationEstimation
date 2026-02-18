function [eta, r] = AnalyticalTPCFwithGSD(model, d, Nr, Rmax, PDF_type, mD, sD)
% Estimates the TPCF for any grain size distribution with given mean mD and
% standard deviation sD of grain diameter for 2D and 3D polycrystals
%
% Inputs:
%  model    : TPCF model, string : "Sha", "Arguelles", "SK"
%  d        : dimension of the medium (2 or 3)
%  Nr       : number of lag distance discretization
%  Rmax     : max lag distance
%  PDF_type : PDF type of grain size distribution
%  mD       : mean grain diameter
%  sD       : standard deviation of grain diameters
%
% Output:
%   eta   :  estimated TPCF
%   r     :  lag distance vector
%
% Example
% [eta, r] = AnalyticalTPCFwithGSD('SK',3,1000,300e-6,'truncated_normal',100e-6,100e-6*0.4);

model = lower(string(model));

r = linspace(0,Rmax,Nr);
eta = zeros(size(r));

[pdf, bounds] = makeGrainDiameterPDF(PDF_type, mD, sD, [1e-6 1-1e-6]);

if d==2

    if model=="sk"
        % Second Moment (Normalization constant for area weighting)
        % E[D^2] = integral( D^2 * p(D) dD )
        ED2 = integral(@(D) D.^2 .* pdf(D), bounds(1), bounds(2));

        % area-Weighted PDF function, p_are(D) = D^2 * p(D) / E[D^2]
        pdf_area = @(D) (D.^2 .* pdf(D)) ./ ED2;
    end

    for i = 1:length(r)
            rr = r(i);
            
            if rr == 0
                eta(i) = 1;
                continue
            end
            
            if rr >= bounds(2)
                eta(i) = 0;
                continue
            end

            if model=="sk"
            
                kernel = @(D) (2/pi)*(acos(rr./D) - (rr./D) .* sqrt( max(1-(rr./D).^2, 0) ));

                integrand = @(D) kernel(D) .* pdf_area(D);
                                
                eta(i) = integral(integrand, max(rr,bounds(1)), bounds(2));
    
            elseif model=="sha"
    
                kernel = @(D) (2/pi)*(acos(rr./D) - (rr./D) .* sqrt( max(1-(rr./D).^2,0) ));

                integrand = @(D) kernel(D) .* pdf(D);
                               
                eta(i) = integral(integrand, max(rr,bounds(1)), bounds(2));
            
            elseif model=="arguelles"
    
                kernel = @(D) exp(-rr./D);
    
                integrand = @(D) kernel(D) .* pdf(D);
    
                eta(i) = integral(integrand, bounds(1), bounds(2));
    
            else
                error('The model should be either Arguelles, Sha, or SK');
            end
    end

elseif d==3

    if model=="sk"
        % Third Moment (Normalization constant for volume weighting)
        % E[D^3] = integral( D^3 * p(D) dD )
        ED3 = integral(@(D) D.^3 .* pdf(D), bounds(1), bounds(2));

        % Volume-Weighted PDF function, p_vol(D) = D^3 * p(D) / E[D^3]
        pdf_vol = @(D) (D.^3 .* pdf(D)) ./ ED3;
    end

    for i = 1:length(r)
        rr = r(i);
        
        if rr == 0
            eta(i) = 1;
            continue
        end
        
        if rr >= bounds(2)
            eta(i) = 0;
            continue;
        end

        if model=="sk"
        
            % Define the Kernel S(r, D) from Eq 14 (or Eq 8 in text)
            % Note: The kernel is 0 if D < r, so we integrate from r to D_max
            kernel = @(D) 1 - 1.5*(rr./D) + 0.5*(rr./D).^3;
            
            integrand = @(D) kernel(D) .* pdf_vol(D);
            
            eta(i) = integral(integrand, max(rr,bounds(1)), bounds(2));

        elseif model=="sha"

            kernel = @(D) 1 - 1.5*(rr./D) + 0.5*(rr./D).^3;

            integrand = @(D) kernel(D) .* pdf(D);

            eta(i) = integral(integrand, max(rr,bounds(1)), bounds(2));
        
        elseif model=="arguelles"

            kernel = @(D) exp(-rr./D);

            integrand = @(D) kernel(D) .* pdf(D);

            eta(i) = integral(integrand, bounds(1), bounds(2));

        else
            error('The model should be either Arguelles, Sha, or SK');
        end

    end

else
    error('Dimension of the problem (d) should be either 2 or 3')
end