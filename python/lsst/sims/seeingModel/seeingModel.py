from builtins import object
from collections import OrderedDict
import numpy as np
from .seeingModelConfig import SeeingModelConfig
from lsst.sims.seeingModel import version


__all__ = ["SeeingModel"]


class SeeingModel(object):
    """LSST FWHM calculations for FWHM_effective and FWHM_geometric.
    Calculations of the delivered values are based on equations in Document-20160
    ("Atmospheric and Delivered Image Quality in OpSim" by Bo Xin, George Angeli, Zeljko Ivezic)
    """
    def configure(self, config=None, configOverrideFile=None):
        """Configure the model. After 'configure' the model config will be frozen.

        Also calculates the fwhm_zenith_system, using self._set_fwhm_zenith_system.

        Parameters
        ----------
        config: SeeingModelConfig, opt
            A configuration class for the seeing model.
            This can can be None, in which case the default values are used.
            The user could override the config parameters before passing,
            or could provide a configOverrideFile.
        configOverrideFile: str, opt
            The filename for a SeeingModelConfig override.
            This file should only contain overrides for the default values.
        """
        if config is None:
            self.config = SeeingModelConfig()
        else:
            if not isinstance(config, SeeingModelConfig):
                raise ValueError('Must use a SeeingModelConfig.')
            self.config = config
        if configOverrideFile is not None:
            self.config.load(configOverrideFile)
        self.config.validate()
        self.config.freeze()
        self._set_fwhm_zenith_system()
        self.filter_list = self.config.filter_list
        self.eff_wavelens = np.array(self.config.filter_effwavelens)

    def efd_requirements(self):
        """Specify data to request from the EFD.

        Returns
        -------
        List of str, float
            List of the EFD columns, delta Time to request information into past (from now).
        """
        return self.config.efd_columns, self.config.efd_delta_time

    def status(self):
        """Report configuration parameters and version information.

        Returns
        -------
        OrderedDict
        """
        status = OrderedDict()
        status['SeeingModel_version'] = '%s' % version.__version__
        status['SeeingModel_sha'] = '%s' % version.__fingerprint__
        for k, v in self.config.iteritems():
            status[k] = v
        return status

    def _set_fwhm_zenith_system(self):
        """Calculate the system contribution to FWHM at zenith.

        This is simply the individual telescope, optics, and camera contributions
        combined in quadrature.
        """
        self.fwhm_system_zenith = np.sqrt(self.config.telescope_seeing**2 +
                                          self.config.optical_design_seeing**2 +
                                          self.config.camera_seeing**2)

    def seeing_at_airmass(self, fwhm_z, airmass=1.0):
        """Calculate the FWHM_eff and FWHM_geom as a function of wavelength, at a range of airmasses,
        given FWHM_z (typically, seeing at 500nm / raw_seeing_wavelength at zenith).

        FWHM_geom represents the geometric size of the PSF; FWHM_eff represents the FWHM of a
        single gaussian which encloses the same number of pixels as N_eff (the number of pixels
        enclosed in the actual PSF -- this is the value to use when calculating SNR).

        FWHM_geom(") = 0.822 * FWHM_eff(") + 0.052"

        The FWHM_eff includes a contribution from the system and from the atmosphere.
        Both of these are expected to scale with airmass^0.6 and with (500(nm)/wavelength(nm))^0.3.
        FWHM_eff = 1.16 * sqrt(FWHM_sys**2 + 1.04*FWHM_atm**2)

        Parameters
        ----------
        fwhm_z : float
            The FWHM_500 (FWHM at 500nm at zenith). Fiducial values is 0.6".
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
        wavelen_correction = np.power(self.config.raw_seeing_wavelength / self.eff_wavelens, 0.3)
        if isinstance(airmass, np.ndarray):
            fwhm_system = self.fwhm_system_zenith * np.outer(np.ones(len(wavelen_correction)),
                                                             airmass_correction)
            fwhm_atmo = fwhm_z * np.outer(wavelen_correction, airmass_correction)
        else:
            fwhm_system = self.fwhm_system_zenith * airmass_correction
            fwhm_atmo = fwhm_z * wavelen_correction * airmass_correction
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
