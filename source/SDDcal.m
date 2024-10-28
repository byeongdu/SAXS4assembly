function SDD = SDDcal(pixel, qvalue, pixelsize, wavelength, whatisinput)
% calculation of SDD from standard sample.
% SDD = SDDcal(pixel, qvalue, pixelsize, wavelength)
% SDD = SDDcal(pixel, qvalue, pixelsize, wavelength, 'pixel')
% px = SDDcal(SDD, qvalue, pixelsize, wavelength, 'SDD')
%   pixel : where the peak of the standard sample occur
%   qvalue : corresponding qvalue for the peak of the standard sample in
%   Angstrom inverse unit.
%   pixelsize : pixelsize of detector in mm unit.
%   wavelength : wavelength in Angstrom unit.
if nargin < 5
    whatisinput = 'pixel';
end
angle = rad2deg(asin(wavelength * qvalue / (4 * pi)));
switch lower(whatisinput)
    case 'sdd'
        SDD = pixel / ( pixelsize / tan(2*deg2rad(angle)));
    case 'pixel'
        SDD = pixel * pixelsize / tan(2*deg2rad(angle));        
end
