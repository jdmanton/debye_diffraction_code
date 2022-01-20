function pupil = get_basic_pupil(sim_params)

px = linspace(-1, 1, sim_params.pupil_size(2));
py = linspace(-1, 1, sim_params.pupil_size(1));
[px, py] = meshgrid(px, py);
r = sqrt(px.^2 + py.^2);
theta = asind(sim_params.numerical_aperture * r / sim_params.refractive_index);


k_0 = 2 * pi / sim_params.wavelength;
k_xy = r * k_0 * sim_params.numerical_aperture;
k_z = sqrt((k_0 * sim_params.refractive_index)^2 - k_xy.^2);


pupil.px = px;
pupil.py = py;
pupil.r = r;
pupil.theta = theta;
pupil.phi = atan2(py, px) * 180 / pi;
pupil.stop = (r <= 1);
pupil.apo = 1 ./ cosd(theta);
pupil.amp = double(pupil.stop);
pupil.phase = double(pupil.stop);
pupil.k_z = k_z;

end
