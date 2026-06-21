function GT = createGroundTruth(vol_size, start_point, end_point, fwhm_axial, fwhm_radial, varargin)

% Input parameters:
% vol_size - Volume size [x_size, y_size, z_size]
% start_point - Start point [x1, y1, z1]
% end_point - End point [x2, y2, z2]
% fwhm_axial - Axial full width at half maximum 
% fwhm_radial - Radial full width at half maximum
%
% Optional parameters (Name-Value pairs):
% 'Intensity' - Peak  intensity (default = 1)
% 'Margin' - Bounding box extension for computation (default = 10)
% 'Precision' - Numerical precision: 'single' or 'double' (default = 'single')
% 'Shape' -  Profile: 'gaussian' or 'uniform' (default = 'gaussian')
%
% Output:
% GT - 3D matrix containing intensity distribution


% =====  Parse input parameters =====
p = inputParser;
addParameter(p, 'Intensity', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Margin', 10, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'Precision', 'single', @(x) ismember(x, {'single', 'double'}));
addParameter(p, 'Shape', 'gaussian', @(x) ismember(x, {'gaussian', 'uniform'}));
parse(p, varargin{:});
params = p.Results;

% ===== Validate inputs =====
validateattributes(vol_size, {'numeric'}, {'positive', 'integer', 'numel', 3});
validateattributes(start_point, {'numeric'}, {'real', 'numel', 3});
validateattributes(end_point, {'numeric'}, {'real', 'numel', 3});
validateattributes(fwhm_axial, {'numeric'}, {'positive', 'scalar'});
validateattributes(fwhm_radial, {'numeric'}, {'positive', 'scalar'});

% ===== Compute  direction =====
direction = end_point - start_point;
length = norm(direction);
if length < eps('double')
    error('Start point and end point cannot be identical.');
end
direction = direction / length;

% ===== Construct local coordinate system =====
% Select a reference vector that is least parallel to the direction
if abs(direction(1)) > abs(direction(2))
    ref_vec = [0, 1, 0];  % Avoid alignment with Y-axis
else
    ref_vec = [1, 0, 0];  % Avoid alignment with X-axis
end

% Compute orthonormal basis vectors
x_dir = cross(ref_vec, direction);
if norm(x_dir) < eps('double')
    x_dir = cross([0, 0, 1], direction); % Backup reference vector
end
x_dir = x_dir / norm(x_dir);
y_dir = cross(direction, x_dir);
y_dir = y_dir / norm(y_dir);

% ===== Compute bounding box (performance optimization) =====
min_point = floor(min(start_point, end_point)) - params.Margin;
max_point = ceil(max(start_point, end_point)) + params.Margin;

% Constrain within volume boundaries
min_point = max(min_point, [1, 1, 1]);
max_point = min(max_point, vol_size);

% Create sub-volume
sub_vol_size = max_point - min_point + 1;
if any(sub_vol_size < 1)
    GT = zeros(vol_size, params.Precision);
    return;
end

% ===== Generate sub-volume grid =====
[X, Y, Z] = meshgrid(...
    min_point(2):max_point(2),...  % Y coordinates
    min_point(1):max_point(1),...  % X coordinates
    min_point(3):max_point(3));    % Z coordinates

if strcmpi(params.Precision, 'single')
    X = single(X);
    Y = single(Y);
    Z = single(Z);
end

% ===== Transform to local coordinates =====
Xc = X - start_point(1);
Yc = Y - start_point(2);
Zc = Z - start_point(3);

x_local = x_dir(1)*Xc + x_dir(2)*Yc + x_dir(3)*Zc;
y_local = y_dir(1)*Xc + y_dir(2)*Yc + y_dir(3)*Zc;
z_local = direction(1)*Xc + direction(2)*Yc + direction(3)*Zc;

% ===== Convert FWHM to Gaussian sigma =====
FWHM_TO_SIGMA = 1/(2*sqrt(2*log(2)));
sigma_axial = fwhm_axial * FWHM_TO_SIGMA;
sigma_radial = fwhm_radial * FWHM_TO_SIGMA;

% Numerical stability protection
sigma_axial = max(sigma_axial, eps(params.Precision));
sigma_radial = max(sigma_radial, eps(params.Precision));

% ===== Compute radial distance =====
r2 = x_local.^2 + y_local.^2;

% ===== Process axial position =====
s = z_local;
s(s < 0) = 0;
s(s > length) = length;

% ===== Generate intensity profile =====
switch params.Shape
    case 'gaussian'
        % Corrected exponential function calculation (multiplication operator added)
        radial_gauss = exp(-r2 / (2 * sigma_radial^2));
        axial_gauss = exp(-(s - length/2).^2 / (2 * sigma_axial^2));
        sub = radial_gauss .* axial_gauss;
        
    case 'uniform'
       % Uniform cylindrical GT with Gaussian axial envelope
        axial_factor = exp(-(s - length/2).^2 / (2 * sigma_axial^2));
        radius_mask = sqrt(r2) <= (fwhm_radial * 0.5); % Radius threshold
        sub = radius_mask .* axial_factor;
end

% ===== Normalize and apply intensity =====
sub = sub / max(sub(:)) * params.Intensity;

% ===== Embed into full volume =====
if strcmpi(params.Precision, 'single')
    GT = zeros(vol_size, 'single');
else
    GT = zeros(vol_size);
end

GT(min_point(1):max_point(1),...
     min_point(2):max_point(2),...
     min_point(3):max_point(3)) = sub;
end