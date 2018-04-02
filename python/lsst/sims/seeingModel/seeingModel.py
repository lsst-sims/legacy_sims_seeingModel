from __future__ import division
from builtins import object
import os
import warnings
import numpy as np
try:
    from lsst.sims.photUtils import Bandpass
    no_photUtils = False
except ImportError:
    no_photUtils = True

__all__ = ["SeeingModel"]


class SeeingModel(object):
    """LSST FWHM calculations for FWHM_effective and FWHM_geometric.
    Calculations of the delivered values are based on equations in Document-20160
    ("Atmospheric and Delivered Image Quality in OpSim" by Bo Xin, George Angeli, Zeljko Ivezic)

    Parameters
    ----------
    telescope_seeing : float, opt
        The contribution to the FWHM at zenith from the telescope, in arcseconds.
        Default 0.25"
    optical_design_seeing : float, opt
        The contribution to the FWHM at zenith from the optical design, in arcseconds.
        Default 0.08"
    camera_seeing : float, opt
        The contribution to the FWHM at zenith from the camera components, in arcseconds.
        Default 0.30"
    raw_seeing_wavelength : float, opt
        The wavelength (in nm) of the provided value of the atmospheric fwhm at zenith.
        Default 500nm.
    filter_effwavelen : numpy.ndarray or None, opt
        An array containing the effective wavelengths per filter in filter_list, in nm.
        If this is None (default), sims_photUtils will be used to calculate the values for ugrizy
        based on the setup throughputs repository.
    """
    def __init__(self, telescope_seeing=0.25,
                 optical_design_seeing=0.08, camera_seeing=0.30,
                 raw_seeing_wavelength=500,
                 filter_effwavelens=None):
        self.set_fwhm_zenith_system(telescope_seeing,
                                    optical_design_seeing,
                                    camera_seeing)
        self.raw_seeing_wavelength = raw_seeing_wavelength
        if filter_effwavelens is None:
            if no_photUtils:

                filter_dict = {'y_effective_wavelength': 971.0,
                               'z_effective_wavelength': 869.1,
                               'u_effective_wavelength': 367.0,
                               'r_effective_wavelength': 622.2,
                               'g_effective_wavelength': 482.5,
                               'i_effective_wavelength': 754.5}

                self.filter_list = ('u', 'g' , 'r', 'i', 'z', 'y')
                self.filter_effwavelens = np.zeros(6, float)

                for i, f in enumerate(self.filter_list):
                    self.filter_effwavelens[i] = filter_dict['%s_effective_wavelength' % f]
                warnings.warn("Could not import sims_photUtils, using defaults", Warning)
            else:
                self._get_effwavelens()
        else:
            self.filter_effwavelens = filter_effwavelens

    def _get_effwavelens(self):
        """Calculate the effective wavelengths, from throughput curves in the $LSST_THOUGHPUTS_DEFAULT dir.
        """
        self.filter_list = ('u', 'g', 'r', 'i', 'z', 'y')
        # Read the throughputs curves from the throughputs package.
        fdir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
        if fdir is None:
            raise ValueError('Set up the throughputs package or define $LSST_THROUGHPUTS_DEFAULT '
                             'to point to the directory containing the throughput curves.')
        lsst = {}
        for f in self.filter_list:
            lsst[f] = Bandpass()
            lsst[f].readThroughput(os.path.join(fdir, 'total_' + f + '.dat'))
        eff_wavelens = np.zeros(len(self.filter_list), float)
        for i, f in enumerate(self.filter_list):
            eff_wavelens[i] = lsst[f].calcEffWavelen()[1]
        self.filter_effwavelens = eff_wavelens

    def set_fwhm_zenith_system(self, telescope_seeing, optical_design_seeing, camera_seeing):
        """Calculate the system contribution to FWHM at zenith.

        This is simply the individual telescope, optics, and camera contributions
        combined in quadrature.
        """
        self.fwhm_system_zenith = np.sqrt(telescope_seeing**2 +
                                          optical_design_seeing**2 +
                                          camera_seeing**2)

    def seeing_at_airmass(self, fwhm_500, airmass=1.0):
        """Calculate the FWHM_eff and FWHM_geom as a function of wavelength, at a range of airmasses,
        given FWHM_500 (seeing at 500nm / raw_seeing_wavelength at zenith).

        FWHM_geom represents the geometric size of the PSF; FWHM_eff represents the FWHM of a
        single gaussian which encloses the same number of pixels as N_eff (the number of pixels
        enclosed in the actual PSF -- this is the value to use when calculating SNR).

        FWHM_geom(") = 0.822 * FWHM_eff(") + 0.052"

        The FWHM_eff includes a contribution from the system and from the atmosphere.
        Both of these are expected to scale with airmass^0.6 and with (500(nm)/wavelength(nm))^0.3.
        FWHM_eff = 1.16 * sqrt(FWHM_sys**2 + 1.04*FWHM_atm**2)

        Parameters
        ----------
        fwhm_500 : float
            The FWHM_500 (FWHM at 500nm at zenith). Fiducial values is 0.6".
        eff_wavelen : numpy.ndarray
            The effective wavelengths of the system bandpasses.
            Can be calculated using get_effwavelens.
        airmass : float or numpy.ndarray
            The airmass at which to calculate the FWHMeff and FWHMgeom values.
            Default 1.0.
            Can be a single value or a numpy array.

        Returns
        -------
        numpy.ndarray, numpy.ndarray
            FWHMeff, FWHMgeom: both are the same shape numpy.ndarray.
            If airmass is a single value, FWHMeff & FWHMgeom are 1-d arrays,
            with the same order as eff_wavelen (i.e. eff_wavelen[0] = u, then FWHMeff[0] = u).
            If airmass is a numpy array, FWHMeff and FWHMgeom are 2-d arrays,
            in the order of <filter><airmass> (i.e. eff_wavelen[0] = u, 1-d array over airmass range).
        """
        airmass_correction = np.power(airmass, 0.6)
        wavelen_correction = np.power(self.raw_seeing_wavelength / self.filter_effwavelens, 0.3)
        if isinstance(airmass, np.ndarray):
            fwhm_system = self.fwhm_system_zenith * np.outer(np.ones(len(wavelen_correction)),
                                                             airmass_correction)
            fwhm_atmo = fwhm_500 * np.outer(wavelen_correction, airmass_correction)
        else:
            fwhm_system = self.fwhm_system_zenith * airmass_correction
            fwhm_atmo = fwhm_500 * wavelen_correction * airmass_correction
        # Calculate combined FWHMeff.
        fwhm_eff = 1.16 * np.sqrt(fwhm_system ** 2 + 1.04 * fwhm_atmo ** 2)
        # Translate to FWHMgeom.
        fwhm_geom = self.fwhmEff_to_fwhmGeom(fwhm_eff)
        return fwhm_eff, fwhm_geom

    @staticmethod
    def fwhmEff_to_fwhmGeom(fwhm_eff):
        """Calculate FWHM_geom from FWHM_eff.

        Parameters
        ----------
        fwhm_eff : float or np.ndarray

        Returns
        -------
        float or np.ndarray
        """
        return (0.822 * fwhm_eff + 0.052)

    @staticmethod
    def fwhmGeom_to_fwhmEff(fwhm_geom):
        """Calculate FWHM_eff from FWHM_geom.

        Parameters
        ----------
        fwhm_geom : float or np.ndarray

        Returns
        -------
        float or np.ndarray
        """
        return (fwhm_geom - 0.052)/0.822
