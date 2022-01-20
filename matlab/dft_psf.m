%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% James Manton, 2020        %
% jmanton@mrc-lmb.cam.ac.uk %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SHOW_IMAGES = 1;
SAVE_IMAGES = 1;
ROOT_NAME = 'dft_psf';

% Simulation parameters
sim_params.wavelength = 500E-9;
sim_params.numerical_aperture = 0.7;
sim_params.refractive_index = 1.33;
sim_params.pupil_size = [256, 256];
sim_params.psf_size = [512, 512, 512] * 1;
sim_params.psf_pitch = [100E-9, 100E-9, 100E-9];

% Set up pupil
pupil = vdc.get_bessel_pupil(sim_params, 0.7, 0.65);
pupil = vdc.apply_polarisation(pupil, 'horizontal');

% Calculate and save PSF
% [electric_field, intensity] = vdc.propagate3d(pupil, sim_params, true);
% if SAVE_IMAGES
%     vdc.save_intensity_16bit(intensity, [ROOT_NAME, '.tif']);
% end

[electric_field, intensity] = vdc.propagate(pupil, 0, sim_params);

figure(1)
imshow(intensity, [])
