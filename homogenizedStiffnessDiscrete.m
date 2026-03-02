function C_eff = homogenizedStiffnessDiscrete(C_rotated, vol_list, material, homoType)
% Computes effective stiffness using Voigt, Reuss, or self-consistent 
% averaging. The latter is only limited to cubic materials.
%
% Inputs
%   C_rotated : 6x6xNg rotated grain stiffness matrices (Pa)
%   vol_list  : Ngx1 grain volumes
%   material  : material struct (elasticConstants in GPa)
%   homoType  : 'voigt', 'reuss', or 'self-consistent'
%
% Output
%   C_eff     : 6x6 effective stiffness matrix (Pa)

    vol = sum(vol_list);
    w = reshape(vol_list, 1, 1, []);
    I6 = eye(6);

    % Always useful, and reused by the self-consistent branch
    C_voigt = (1/vol) * squeeze(sum(C_rotated .* w, 3));

    switch lower(char(homoType))

        case {'voigt','v'}
            C_eff = C_voigt;

        case {'reuss','r'}
            Ng = size(C_rotated, 3);
            S_rotated = zeros(6,6,Ng);

            for ig = 1:Ng
                % More stable than inv(C_rotated(:,:,ig))
                S_rotated(:,:,ig) = C_rotated(:,:,ig) \ I6;
            end

            S_eff = (1/vol) * squeeze(sum(S_rotated .* w, 3));
            C_eff = S_eff \ I6;

        case {'sc','selfconsistent','self-consistent'}
            % Implemented here for cubic crystals only, following Norouzian 2019a
            if numel(material.elasticConstants) ~= 3
                error('Self-consistent averaging is implemented here only for cubic crystals.');
            end

            c11 = material.elasticConstants(1); % GPa
            c12 = material.elasticConstants(2); % GPa
            c44 = material.elasticConstants(3); % GPa

            nu_c = c11 - c12 - 2*c44; % second-order anisotropy coefficient (GPa)

            % If the crystal is effectively isotropic, all schemes collapse
            tol = 1e-12 * max(abs([c11, c12, c44, 1]));
            if abs(nu_c) < tol
                disp('The given material is isotropic: all schemes collapse')
                C_eff = C_voigt;
            else
                % Isotropic part of the local cubic stiffness: C^I
                % Note : c11^iso = c12 + 2*c44;  
                C_iso = stiffnessMatrix('iso', c12 + 2*c44, c44) * 1e9; % Pa

                % Solve for h from Eqs. (15)-(16)
                h = selfConsistentCubicParam(c11, c12, c44);

                % Isotropic part C^{I,sc} from Eq. (19)
                lambda_sc = h*c11 + (1-h)*c12;   % GPa
                mu_sc     = (1 + 2*h)*c44;       % GPa
                C_iso_sc  = stiffnessMatrix('iso', lambda_sc + 2*mu_sc, mu_sc) * 1e9; % Pa

                % Eq. (23c): C_eff^sc = C^{I,sc} + c_sc * <R>
                % Since C_voigt = C^I + nu_c * <R>, we use:
                % <R> = (C_voigt - C_iso_local) / nu_c
                c_sc = nu_c*(1 - 3*h) - 10*h*c44; % GPa

                C_eff = C_iso_sc + (c_sc / nu_c) * (C_voigt - C_iso);
            end

        otherwise
            error('Unknown homoType "%s". Use ''voigt'', ''reuss'', or ''self-consistent''.', char(homoType));
    end

    % Enforce symmetry numerically
    C_eff = 0.5 * (C_eff + C_eff.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Helper function for self consistent technique
function h = selfConsistentCubicParam(c11, c12, c44)
% Compute c44^sc and h for the cubic self-consistent scheme.
%
% Inputs are in GPa. Output c44_sc is in GPa, h is dimensionless.

    % Eq. (16):
    % 8*(c44_sc)^3 + (5*c11 + 4*c12)*(c44_sc)^2
    % - c44*(7*c11 - 4*c12)*c44_sc
    % - c44*(c11 - c12)*(c11 + 2*c12) = 0
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
end