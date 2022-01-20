import vdc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sim_params = {
    'pupil_size': np.array([256, 256]),
    'psf_size': np.array([256, 256, 256]),
    'numerical_aperture': 1.2,
    'refractive_index': 1.33,
    'wavelength': 488E-9,
    'psf_pitch': np.array([50E-9, 50E-9, 50E-9])
}

pupil = vdc.Pupil(sim_params)
pupil.set_bessel_pupil(sim_params, 0.9, 0.8)
pupil.apply_polarisation('vertical')
electric_field, intensity = pupil.propagate(0, sim_params)

e3d, i3d = pupil.propagate3d(sim_params)

fig, ax = plt.subplots(2, 2)
ax[0, 0].imshow(pupil.amp)
ax[0, 0].title.set_text('Pupil amplitude')
ax[0, 1].imshow(np.angle(pupil.phase))
ax[0, 1].title.set_text('Pupil phase')
ax[1, 0].imshow(intensity)
ax[1, 0].title.set_text('PSF')
ax[1, 1].imshow(np.abs(np.fft.fftshift(np.fft.fft2(intensity))))
ax[1, 1].title.set_text('MTF')
plt.tight_layout()
plt.show()
