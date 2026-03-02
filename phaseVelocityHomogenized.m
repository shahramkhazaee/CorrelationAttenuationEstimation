function [V, Ceff] = phaseVelocityHomogenized(Cloc, rho, homogType, tol)
% phaseVelocityHomogenized
% Wrapper that combines Voigt, Reuss, and Self-Consistent (cubic only) into one function.
%
% Inputs:
%   Cloc      [6x6] stiffness in Voigt notation (11,22,33,23,13,12), Pa
%   rho       density, kg/m^3
%   homogType 'voigt' | 'reuss' | 'hill' | 'sc'
%   tol       relative tolerance for cubic check (default 1e-6)
%
% Outputs:
%   V    [2x1] Vp, Vs (m/s) for the selected isotropic effective media
%   Ceff [6x6] isotropic homogenized background stiffness matrix
%   info struct with some intermediate values (optional)
%
% References:
% - R. deWit (2008), Elastic constants and thermal expansion averages of a 
%   nontextured polycrystal
%
%  - M. Norouzian & J.A. Turner (2019) "Ultrasonic wave propagation predictions
%    for polycrystalline materials using threedimensional synthetic microstructures:
%    Phase velocity variations" (JASA 145, 2171).

    if nargin < 3 || isempty(homogType)
        homogType = 'voigt';
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-6;
    end

    if ~isequal(size(Cloc), [6 6])
        error('Cloc must be a 6x6 stiffness matrix in Voigt notation.');
    end
    if ~isscalar(rho) || ~isfinite(rho) || rho <= 0
        error('rho must be a positive finite scalar.');
    end
    if ~isscalar(tol) || ~isfinite(tol) || tol <= 0
        error('tol must be a positive finite scalar.');
    end


    switch lower(char(homogType))

        case 'voigt'
            % Voigt (uniform strain) isotropic aggregate velocities from a 
            % general 6x6 stiffness C.

            % Enforce symmetry
            C = 0.5 * (Cloc + Cloc.');

            A = C(1,1) + C(2,2) + C(3,3);
            B = C(1,2) + C(1,3) + C(2,3);
            D = C(4,4) + C(5,5) + C(6,6);

            % Eq. (15) of R. deWit
            KV = (A + 2*B) / 9;       % Voigt bulk modulus
            GV = (A - B + 3*D) / 15;  % Voigt shear modulus 

            % Eq. (6) of R. deWit
            C11iso = KV + (4/3)*GV;
            C44iso = GV;
            Ceff = stiffnessMatrix('iso', C11iso, C44iso);

            if KV <= 0 || GV <= 0
                error('Computed Voigt moduli are non-positive. Check C units/convention.');
            end

            V = [sqrt(C11iso/rho); sqrt(C44iso/rho)];

        case 'reuss'
            % Reuss (uniform stress) isotropic aggregate velocities from a 
            % general 6x6 stiffness C.

            % Enforce symmetry
            C = 0.5 * (Cloc + Cloc.');

            % Compliance S = inv(C)
            S = C \ eye(6);

            A = S(1,1) + S(2,2) + S(3,3);
            B = S(1,2) + S(1,3) + S(2,3);
            D = S(4,4) + S(5,5) + S(6,6);

            % Eq. (18) of R. deWit
            KR = 1 / (A + 2*B);         % Reuss bulk modulus
            GR = 15 / (4*(A-B) + 3*D);  % Reuss shear modulus

            if KR <= 0 || GR <= 0
                error('Computed Reuss moduli are non-positive. Check C units/convention.');
            end

            % Eq. (6) of R. deWit
            C11iso = KR + (4/3)*GR;
            C44iso = GR;
            Ceff = stiffnessMatrix('iso', C11iso, C44iso);

            V = [sqrt(C11iso/rho); sqrt(C44iso/rho)];

        case 'hill'
            % Average between Voigt & Reuss averages

            % Enforce symmetry
            C = 0.5 * (Cloc + Cloc.');

            % Voigt average
            A = C(1,1) + C(2,2) + C(3,3);
            B = C(1,2) + C(1,3) + C(2,3);
            D = C(4,4) + C(5,5) + C(6,6);

            % Eq. (15) of R. deWit
            KV = (A + 2*B) / 9;       % Voigt bulk modulus
            GV = (A - B + 3*D) / 15;  % Voigt shear modulus 

            if KV <= 0 || GV <= 0
                error('Computed Voigt moduli are non-positive. Check C units/convention.');
            end

            CeffV = stiffnessMatrix('iso', KV + (4/3)*GV, GV);
            
            % Reuss average
            S = C \ eye(6);

            A = S(1,1) + S(2,2) + S(3,3);
            B = S(1,2) + S(1,3) + S(2,3);
            D = S(4,4) + S(5,5) + S(6,6);

            % Eq. (18) of R. deWit
            KR = 1 / (A + 2*B);         % Reuss bulk modulus
            GR = 15 / (4*(A-B) + 3*D);  % Reuss shear modulus

            if KR <= 0 || GR <= 0
                error('Computed Reuss moduli are non-positive. Check C units/convention.');
            end

            CeffR = stiffnessMatrix('iso', KR + (4/3)*GR, GR);

            Ceff = 0.5 * (CeffV + CeffR);

            V = [sqrt(Ceff(1,1)/rho); sqrt(Ceff(4,4)/rho)];

        case {'sc','selfconsistent','self-consistent'}
            % Self-consistent isotropic equivalent for a random polycrystal of cubic grains
            % using the notation and equations of Norouzian & Turner (JASA 2019).

            % Symmetrize for numerical cleanliness
            C = 0.5 * (Cloc + Cloc.');

            ok = isCubic(C, tol);
            if ~ok
                error('Cloc is not cubic (in the standard cubic form).');
            end

            % Extract cubic constants
            c11 = C(1,1);
            c12 = C(1,2);
            c44 = C(4,4);

            % Solve for c44_sc from Eq. (16) (lone positive root)
            % Eq. (16): 8 x^3 + (5c11+4c12) x^2 - c44(7c11-4c12) x - c44(c11-c12)(c11+2c12) = 0
            p = [ 8, ...
                  (5*c11 + 4*c12), ...
                  -c44*(7*c11 - 4*c12), ...
                  -c44*(c11 - c12)*(c11 + 2*c12) ];
            
            rts = roots(p);

            % Keep the only positive root
            c44_sc = max(real(rts));

            % Compute h from Eq. (15) of Nourouzian & Turner 2019
            % Eq. (15): h = ((c11+2c12+6c44_sc)(c11-c12-2c44_sc)) / 
            % (3*(8c44_sc^2 + 9c11 c44_sc + (c11-c12)(c11+2c12)))
            num = (c11 + 2*c12 + 6*c44_sc) * (c11 - c12 - 2*c44_sc);
            den = 3 * ( 8*(c44_sc^2) + 9*c11*c44_sc + (c11 - c12)*(c11 + 2*c12) );
            h = num / den;

            % Compute isotropic infinite-polycrystal SC average from Eq.
            % (21) of Nourouzian & Turner 2019
            
            % Second-order anisotropy coefficient (Eq. (2)): 
            nu = c11 - c12 - 2*c44;

            % Eq. (21) of Nourouzian & Turner 2019 gives the isotropic 
            % Lamé-type coefficients directly:
            % lambda_sc = C12_sc = c12 + (c/5)*(1 + 2h)
            % mu_sc     = C44_sc = c44 + (c/5)*(1 - 3h)
            C12_sc = c12 + (nu/5) * (1 + 2*h);
            C44_sc = c44 + (nu/5) * (1 - 3*h);
            C11_sc = C12_sc + 2*C44_sc; % isotropy identity

            if C44_sc <= 0 || C11_sc <= 0
                error('SC produced non-positive isotropic constants. Check Cloc stability.');
            end

            % Build isotropic stiffness matrix Ciso (Voigt)
            Ceff = stiffnessMatrix('iso',C11_sc, C44_sc);

            V = [sqrt(C11_sc/rho); sqrt(C44_sc/rho)];

        otherwise
            error('Unknown homogType. Use ''voigt'', ''reuss'', or ''sc''.');
    end

end

function tf = isCubic(C, tol)
% isCubic  Check if a 6x6 stiffness matrix is cubic in the standard cubic axes.
%
% Voigt order assumed: (11,22,33,23,13,12).
%
% tf = true if cubic conditions hold within tolerance.
% tol is a relative tolerance (default 1e-6).

    if nargin < 2 || isempty(tol)
        tol = 1e-6;
    end

    if ~isequal(size(C), [6 6])
        error('C must be 6x6.');
    end
    if ~isscalar(tol) || tol <= 0 || ~isfinite(tol)
        error('tol must be a positive finite scalar.');
    end

    % Scale for relative comparisons
    scale = max(1, max(abs(C(:))));
    rel = @(x) abs(x) / scale;

    % Symmetry
    symErr = rel(C - C.');

    % Groups that must be equal
    c11s = [C(1,1), C(2,2), C(3,3)];
    c12s = [C(1,2), C(1,3), C(2,3)];
    c44s = [C(4,4), C(5,5), C(6,6)];

    err_c11 = max(abs(c11s - mean(c11s))) / scale;
    err_c12 = max(abs(c12s - mean(c12s))) / scale;
    err_c44 = max(abs(c44s - mean(c44s))) / scale;

    % Entries that must be ~0 in standard cubic form
    normalShear = C(1:3, 4:6);
    shearShearOff = [C(4,5), C(4,6), C(5,6)];

    err_normalShear = max(abs(normalShear(:))) / scale;
    err_shearShearOff = max(abs(shearShearOff(:))) / scale;

    % Combine decisions
    tf = (max(symErr(:)) <= tol) && ...
         (err_c11 <= tol) && (err_c12 <= tol) && (err_c44 <= tol) && ...
         (err_normalShear <= tol) && (err_shearShearOff <= tol);
end