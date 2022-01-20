function pupil = get_bessel_pupil(sim_params, outer_na, inner_na)

outer_r = outer_na / sim_params.numerical_aperture;
inner_r = inner_na / sim_params.numerical_aperture;
pupil = vdc.get_basic_pupil(sim_params);
pupil.amp = (pupil.r <= outer_r & pupil.r >= inner_r);

end