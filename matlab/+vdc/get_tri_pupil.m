function pupil = get_tri_pupil(sim_params, outer_na, inner_na)

outer_r = outer_na / sim_params.numerical_aperture;
inner_r = inner_na / sim_params.numerical_aperture;
pupil = vdc.get_basic_pupil(sim_params);
pupil.amp = (pupil.r <= outer_r & pupil.r >= inner_r);

band_radius = 0.5 * (outer_r + inner_r);
band_x = 0.5 * band_radius .* cosd(60);
band_mask = (round(sim_params.psf_size(1) * pupil.px .* 0.5) == round(sim_params.psf_size(1) * band_x));
band_mask = band_mask + (round(sim_params.psf_size(1) * pupil.px .* 0.5) == round(-sim_params.psf_size(1) * band_radius .* 0.5));
pupil.amp = pupil.amp .* band_mask;

end