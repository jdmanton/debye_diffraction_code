function pupil = apply_polarisation(pupil, type)

% See 'Inversion of the Debye-Wolf diffraction integral using an
% eigenfunction representation of the electric fields in the focal region'
% by M. R. Foreman, S. S. Sherif, P. R. T. Munro, and P. Török 
% https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-16-7-4901&id=156860

% L = [(1 + cosd(pupil.theta)) - (1 - cosd(pupil.theta) .* cosd(2 .* pupil.phi)), ...
%          -(1 - cosd(pupil.theta)) .* sind(2 .* pupil.phi); ...
%      -(1 - cosd(pupil.theta)) .* sind(2 * pupil.phi), ...
%          (1 + cosd(pupil.theta)) + (1 - cosd(pupil.theta)) .* cosd(2 * pupil.phi); ...
%      -2 * sind(pupil.theta) .* cosd(pupil.phi), ...
%          -2 * sind(pupil.theta) .* cosd(pupil.phi)];
%
% e = L * [pupil.ex, pupil.ey];

switch type
    case 'vertical'
%         pupil.pol1 = ((1 + cosd(pupil.theta)) - (1 - cosd(pupil.theta)) .* cosd(2 * pupil.phi));
%         pupil.pol2 = -(1 - cosd(pupil.theta)) .* sind(2 * pupil.phi);
%         pupil.pol3 = -2 * cosd(pupil.phi) .* sind(pupil.theta);
        x = 0;
        y = 1;
    case 'horizontal'
        x = 1;
        y = 0;
    case 'circular'
        x = 1;
        y = 1i;
    case 'radial'
        x = cosd(pupil.phi);
        y = sind(pupil.phi);
    case 'azimuthal'
        x = -sind(pupil.phi);
        y = cosd(pupil.phi);
    case 'dipole_x'
%         x = (cosd(pupil.theta) - 1) .* sind(pupil.phi) .* cosd(pupil.phi);
%         y = cosd(pupil.theta) .* sind(pupil.phi).^2 + cosd(pupil.phi).^2;
        
        x = cosd(pupil.theta) .* cosd(pupil.phi).^2 + sind(pupil.phi).^2;
        y = (cosd(pupil.theta) - 1) .* sind(pupil.phi) .* cosd(pupil.phi);
    case 'dipole_y'
%         x = cosd(pupil.theta) .* cosd(pupil.phi).^2 + sind(pupil.phi).^2;
%         y = (cosd(pupil.theta) - 1) .* sind(pupil.phi) .* cosd(pupil.phi);

        x = (cosd(pupil.theta) - 1) .* sind(pupil.phi) .* cosd(pupil.phi);
        y = cosd(pupil.theta) .* sind(pupil.phi).^2 + cosd(pupil.phi).^2;
    case 'dipole_z'
        x = sind(pupil.theta) .* cosd(pupil.phi);
        y = sind(pupil.theta) .* sind(pupil.phi);
    otherwise
        error(['Polarisation type ''', type, ''' not implemented'])
end

pupil.pol1 = x .* (1 - cosd(2 .* pupil.phi) .* (1 - cosd(pupil.theta)) + cosd(pupil.theta)) + y .* (-1 + cosd(pupil.theta)) .* sind(2 .* pupil.phi);
pupil.pol2 = y .* (1 + cosd(2 .* pupil.phi) .* (1 - cosd(pupil.theta)) + cosd(pupil.theta)) + x .* (-1 + cosd(pupil.theta)) .* sind(2 .* pupil.phi);
pupil.pol3 = -2 .* x .* cosd(pupil.phi) .* sind(pupil.theta) - 2 .* y .* sind(pupil.phi) .* sind(pupil.theta);

pupil.ex = pupil.pol1 .* pupil.stop .* pupil.apo .* pupil.amp .* pupil.phase;
pupil.ey = pupil.pol2 .* pupil.stop .* pupil.apo .* pupil.amp .* pupil.phase;
pupil.ez = pupil.pol3 .* pupil.stop .* pupil.apo .* pupil.amp .* pupil.phase;

end
