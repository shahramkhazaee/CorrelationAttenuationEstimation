function C = stiffnessMatrix(symmetry, varargin)
% Build a 6x6 stiffness matrix (Voigt notation: 11,22,33,23,13,12)
%
% Supported symmetries:
%   iso:        C = stiffnessMatrix('iso', c11, c44)
%   cubic:      C = stiffnessMatrix('cubic', c11, c12, c44)
%   hex/ti:     C = stiffnessMatrix('hex', c11, c12, c13, c33, c44)
%   orthotropic:C = stiffnessMatrix('ortho', c11,c22,c33,c12,c13,c23,c44,c55,c66)
%   tetragonal: C = stiffnessMatrix('tetragonal', c11,c12,c13,c33,c44,c66)
%   trigonal:   C = stiffnessMatrix('trigonal', c11,c12,c13,c14,c33,c44)
%
% Notes:
% - iso expects (c11,c44) and sets c12 = c11 - 2*c44.
% - hex and ti are equivalent in this implementation.

symmetry = lower(strtrim(char(symmetry)));

if nargin < 2
    error('At least two arguments required: symmetry type and stiffness parameters.');
end

% Support for symbolic or numeric
if isa(varargin{1},'sym')
    C = sym(zeros(6,6));
else
    C = zeros(6,6);
end

switch symmetry

    case {'iso','isotropic'}
        if numel(varargin) ~= 2
            error('Isotropic: (c11, c44) required.');
        end
        c11 = varargin{1};
        c44 = varargin{2};
        c12 = c11 - 2*c44;

        C(1,1)=c11; C(2,2)=c11; C(3,3)=c11;
        C(1,2)=c12; C(1,3)=c12; C(2,3)=c12;
        C(2,1)=c12; C(3,1)=c12; C(3,2)=c12;
        C(4,4)=c44; C(5,5)=c44; C(6,6)=c44;

    case 'cubic'
        if numel(varargin) ~= 3
            error('Cubic: (c11, c12, c44) required.');
        end
        c11 = varargin{1};
        c12 = varargin{2};
        c44 = varargin{3};

        C(1,1)=c11; C(2,2)=c11; C(3,3)=c11;
        C(1,2)=c12; C(1,3)=c12; C(2,3)=c12;
        C(2,1)=c12; C(3,1)=c12; C(3,2)=c12;
        C(4,4)=c44; C(5,5)=c44; C(6,6)=c44;

    case {'hex','hexagonal','ti','transverselyisotropic','transversely-isotropic'}

        if numel(varargin) ~= 5
            error('Hexagonal: (c11, c12, c13, c33, c44) required.');
        end

        c11 = varargin{1}; c12 = varargin{2}; c13 = varargin{3};
        c33 = varargin{4}; c44 = varargin{5};
        c66 = (c11 - c12)/2;

        C(1,1) = c11; C(2,2) = c11; C(3,3) = c33;
        C(1,2) = c12; C(2,1) = c12;
        C(1,3) = c13; C(3,1) = c13;
        C(2,3) = c13; C(3,2) = c13;
        C(4,4) = c44; C(5,5) = c44;
        C(6,6) = c66;

    case {'ortho','orthotropic'}
        if numel(varargin) ~= 9
            error('Orthotropic: (c11,c22,c33,c12,c13,c23,c44,c55,c66) required.');
        end
        c11=varargin{1}; c22=varargin{2}; c33=varargin{3};
        c12=varargin{4}; c13=varargin{5}; c23=varargin{6};
        c44=varargin{7}; c55=varargin{8}; c66=varargin{9};

        C(1,1)=c11; C(2,2)=c22; C(3,3)=c33;
        C(1,2)=c12; C(2,1)=c12;
        C(1,3)=c13; C(3,1)=c13;
        C(2,3)=c23; C(3,2)=c23;
        C(4,4)=c44; C(5,5)=c55; C(6,6)=c66;

    case 'tetragonal'
        if numel(varargin) ~= 6
            error('Tetragonal: (c11,c12,c13,c33,c44,c66) required.');
        end
        c11=varargin{1}; c12=varargin{2}; c13=varargin{3};
        c33=varargin{4}; c44=varargin{5}; c66=varargin{6};

        C(1,1)=c11; C(2,2)=c11; C(3,3)=c33;
        C(1,2)=c12; C(2,1)=c12;
        C(1,3)=c13; C(3,1)=c13;
        C(2,3)=c13; C(3,2)=c13;
        C(4,4)=c44; C(5,5)=c44; C(6,6)=c66;

    case 'trigonal'
        if numel(varargin) ~= 6
            error('Trigonal: (c11,c12,c13,c14,c33,c44) required.');
        end
        c11=varargin{1}; c12=varargin{2}; c13=varargin{3};
        c14=varargin{4}; c33=varargin{5}; c44=varargin{6};
        c66 = (c11 - c12)/2;

        % Standard trigonal form (class 32) in the standard axes
        C(1,1)=c11; C(2,2)=c11; C(3,3)=c33;
        C(1,2)=c12; C(2,1)=c12;
        C(1,3)=c13; C(3,1)=c13;
        C(2,3)=c13; C(3,2)=c13;

        C(1,4)=c14;  C(4,1)=c14;
        C(2,4)=-c14; C(4,2)=-c14;
        C(5,6)=c14;  C(6,5)=c14;

        C(4,4)=c44; C(5,5)=c44; C(6,6)=c66;

    otherwise
        error('Unknown symmetry type: use ''iso'', ''cubic'', ''hex'', ''ortho'', ''tetragonal'', or ''trigonal''.');
end
