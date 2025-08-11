import numpy as np
from tqdm import tqdm

def dft2(X, k, a, N):
    f1 = np.linspace(-1 / (2 * a[0]), 1 / (2 * a[0]), N[0]) - k[0] / X.shape[0]
    f1 = f1.reshape((-1, 1))
    f2 = np.linspace(-1 / (2 * a[1]), 1 / (2 * a[1]), N[1]) - k[1] / X.shape[1]
    x1 = np.arange(0, X.shape[0])
    x2 = np.arange(0, X.shape[1])
    x2 = x2.reshape((-1, 1))
    F1 = np.exp(-1j * 2 * np.pi * f1 * x1)
    F2 = np.exp(-1j * 2 * np.pi * x2 * f2)
    Xhat = np.matmul(F1, np.matmul(X, F2))
    return Xhat


class Pupil:
    def __init__(self, sim_params):
        px = np.linspace(-1, 1, sim_params['pupil_size'][0])
        py = np.linspace(-1, 1, sim_params['pupil_size'][1])
        self.px, self.py = np.meshgrid(px, py)
        self.r = np.sqrt(self.px**2 + self.py**2)
        self.r[self.r > 1] = 0
        self.theta = np.arcsin(sim_params['numerical_aperture'] * self.r / sim_params['refractive_index'])
        self.phi = np.arctan2(self.py, self.px) * 180 / np.pi
        self.stop = (np.sqrt(self.px**2 + self.py**2) <= 1)
        self.apo = 1 / np.cos(self.theta)
        self.amp = self.stop.astype(float)
        self.phase = self.stop.astype(float)
        self.k_0 = 2 * np.pi / sim_params['wavelength']
        self.k_xy = self.r * self.k_0 * sim_params['numerical_aperture']
        self.k_z = np.sqrt((self.k_0 * sim_params['refractive_index'])**2 - self.k_xy**2)


    def set_bessel_pupil(self, sim_params, outer_na, inner_na):
        outer_r = outer_na / sim_params['numerical_aperture']
        inner_r = inner_na / sim_params['numerical_aperture']
        self.amp = (self.r <= outer_r) & (self.r >= inner_r)
    

    def propagate(self, z, sim_params):
        scaling = sim_params['wavelength'] / (2 * sim_params['numerical_aperture'] * sim_params['psf_pitch'][0:2])
        defocus = np.exp(1j * self.k_z * z)

        field_x = dft2(defocus * self.ex, np.array([0, 0]), scaling, sim_params['psf_size'][0:2])
        field_y = dft2(defocus * self.ey, np.array([0, 0]), scaling, sim_params['psf_size'][0:2])
        field_z = dft2(defocus * self.ez, np.array([0, 0]), scaling, sim_params['psf_size'][0:2])

        electric_field = np.zeros([sim_params['psf_size'][0], sim_params['psf_size'][1], 3], dtype=complex)

        electric_field[:, :, 0] = field_x
        electric_field[:, :, 1] = field_y
        electric_field[:, :, 2] = field_z

        intensity = np.abs(field_x)**2 + np.abs(field_y)**2 + np.abs(field_z)**2

        return electric_field, intensity
    
    
    def propagate3d(self, sim_params):
        zpx = np.linspace(-(sim_params['psf_size'][2] - 1) / 2, (sim_params['psf_size'][2] - 1) / 2, sim_params['psf_size'][2])
        electric_field = np.zeros([sim_params['psf_size'][0], sim_params['psf_size'][1], sim_params['psf_size'][2], 3], dtype=complex)
        intensity = np.zeros(sim_params['psf_size'])

        for z_index in tqdm(range(1, sim_params['psf_size'][2])):
            z = zpx[z_index] * sim_params['psf_pitch'][2]
            electric_field[:, :, z_index, :], intensity[:, :, z_index] = self.propagate(z, sim_params)
        
        return electric_field, intensity
    

    def apply_polarisation(self, polarisation):
        x, y = 0, 0
        if (polarisation == 'vertical'):
            x, y = 0, 1
        elif (polarisation == 'horizontal'):
            x, y = 1, 0
        elif (polarisation == 'circular'):
            x, y = 1, 1j
        elif(polarisation == 'radial'):
            x, y = np.cos(self.phi), np.sin(self.phi)
        elif (polarisation == 'azimuthal'):
            x, y = -np.sin(self.phi), np.cos(self.phi)
        elif (polarisation == 'dipole_x'):
            x = np.cos(self.theta) * np.cos(self.phi)**2 + np.sin(self.phi)**2
            y = (np.cos(self.theta) - 1) * np.sin(self.phi) * np.cos(self.pupil.phi)
        elif (polarisation == 'dipole_y'):
            x = (np.cos(self.theta) - 1) * np.sin(self.phi) * np.cos(self.pupil.phi)
            y = np.cos(self.theta) * np.cos(self.phi)**2 + np.sin(self.phi)**2
        elif (polarisation == 'dipole_z'):
            x = np.sin(self.theta) * np.cos(self.phi)
            y = np.sin(self.theta) * np.sin(self.phi)
        
        self.pol1 = x * (1 - np.cos(2 * self.phi) * (1 - np.cos(self.theta)) + np.cos(self.theta)) + y * (-1 + np.cos(self.theta)) * np.sin(2 * self.phi)
        self.pol2 = y * (1 - np.cos(2 * self.phi) * (1 - np.cos(self.theta)) + np.cos(self.theta)) + x * (-1 + np.cos(self.theta)) * np.sin(2 * self.phi)
        self.pol3 = -2 * x * np.cos(self.phi) * np.sin(self.theta) - 2 * y * np.sin(self.phi) * np.sin(self.theta)

        self.ex = self.pol1 * self.stop * self.apo * self.amp * np.exp(1j * self.phase)
        self.ey = self.pol2 * self.stop * self.apo * self.amp * np.exp(1j * self.phase)
        self.ez = self.pol3 * self.stop * self.apo * self.amp * np.exp(1j * self.phase)
