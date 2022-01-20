function [electric_field, intensity] = propagate3d(pupil, sim_params, varargin)

if nargin > 2
    DISPLAY_PROGRESS = varargin{1};
else
    DISPLAY_PROGRESS = 0;
end

zpx = linspace(-(sim_params.psf_size(3) - 1) / 2, (sim_params.psf_size(3) - 1) / 2, sim_params.psf_size(3));
electric_field = zeros(sim_params.psf_size(1), sim_params.psf_size(2), sim_params.psf_size(3), 3);
intensity = zeros(sim_params.psf_size(1), sim_params.psf_size(2), sim_params.psf_size(3));

if DISPLAY_PROGRESS
    progress = waitbar(0, 'Calculating focus...');
end

for z_index = 1:sim_params.psf_size(3)
    z = zpx(z_index) * sim_params.psf_pitch(3);
    [electric_field(:, :, z_index, :), intensity(:, :, z_index)] = vdc.propagate(pupil, z, sim_params);
    if DISPLAY_PROGRESS
        waitbar(z_index / sim_params.psf_size(3))
    end
end

if DISPLAY_PROGRESS
    close(progress)
end

end
