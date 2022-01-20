function [electric_field, intensity] = propagate(pupil, z, sim_params)

scaling = sim_params.wavelength ./ (2 * sim_params.numerical_aperture * sim_params.psf_pitch);

defocus = exp(1i * pupil.k_z * z);

field_x = vdc.dft2(defocus .* pupil.ex, [0, 0], scaling, sim_params.psf_size(1:2));
field_y = vdc.dft2(defocus .* pupil.ey, [0, 0], scaling, sim_params.psf_size(1:2));
field_z = vdc.dft2(defocus .* pupil.ez, [0, 0], scaling, sim_params.psf_size(1:2));

electric_field(:, :, 1) = field_x;
electric_field(:, :, 2) = field_y;
electric_field(:, :, 3) = field_z;

intensity = abs(field_x).^2 + abs(field_y).^2 + abs(field_z).^2;

end
