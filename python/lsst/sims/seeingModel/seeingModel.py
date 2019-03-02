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

    Parameters
    ----------
    config: SeeingModelConfig, opt
        A configuration class for the seeing model.
        This can be None, in which case the default SeeingModelConfig is used.
        The user should set any non-default values for SeeingModelConfig before
        configuration of the actual SeeingModel.

    self.efd_requirements and self.map_requirements are also set.
    efd_requirements is a tuple: (list of str, float).
    This corresponds to the data columns required from the EFD and the amount of time history required.
    map_requirements is a list of str.
    This corresponds to the data columns required in the map dictionary passed when calculating the
    processed telemetry values.
    """
    def __init__ (self, config=None):
        self._configure(config=config)
        self.efd_requirements = (self._config.efd_columns, self._config.efd_delta_time)
        self.map_requirements = ['airmass']
        self.seeingcol = self._config.efd_columns[0]

    def _configure(self, config=None):
        """Configure the model. After 'configure' the model config will be frozen.

        Also calculates the fwhm_zenith_system, using self._set_fwhm_zenith_system.

        Parameters
        ----------
        config: SeeingModelConfig, opt
            A configuration class for the seeing model.
            This can be None, in which case the default values are used.
        """
        if config is None:
            self._config = SeeingModelConfig()
        else:
            if not isinstance(config, SeeingModelConfig):
                raise ValueError('Must use a SeeingModelConfig.')
            self._config = config
        self._config.validate()
        self._config.freeze()
        self._set_fwhm_zenith_system()
        self.filter_list = tuple(self._config.filter_list)
        self.eff_wavelens = np.array(self._config.filter_effwavelens)

    def status(self):
        """Report configuration parameters and version information.

        Returns
        -------
        OrderedDict
        """
        status = OrderedDict()
        status['SeeingModel_version'] = '%s' % version.__version__
        status['SeeingModel_sha'] = '%s' % version.__fingerprint__
        for k, v in self._config.iteritems():
            status[k] = v
        status['map_columns'] = self.map_requirements
        return status

    def _set_fwhm_zenith_system(self):
        """Calculate the system contribution to FWHM at zenith.

        This is simply the individual telescope, optics, and camera contributions
        combined in quadrature.
        """
        self.fwhm_system_zenith = np.sqrt(self._config.telescope_seeing**2 +
                                          self._config.optical_design_seeing**2 +
                                          self._config.camera_seeing**2)

    def __call__(self, efdData, mapDict):
        """Calculate the seeing values FWHM_eff and FWHM_geom at the given airmasses,
        for the specified effective wavelengths, given FWHM_zenith (typically FWHM_500).

        FWHM_geom represents the geometric size of the PSF; FWHM_eff represents the FWHM of a
        single gaussian which encloses the same number of pixels as N_eff (the number of pixels
        enclosed in the actual PSF -- this is the value to use when calculating SNR).

        FWHM_geom(") = 0.822 * FWHM_eff(") + 0.052"

        The FWHM_eff includes a contribution from the system and from the atmosphere.
        Both of these are expected to scale with airmass^0.6 and with (500(nm)/wavelength(nm))^0.3.
        FWHM_eff = 1.16 * sqrt(FWHM_sys**2 + 1.04*FWHM_atm**2)

        Parameters
        ----------
        efdData: dict
            Dictionary of input telemetry, typically from the EFD.
            This must contain columns self.efd_requirements.
            (work in progress on handling time history).
        mapDict: dict
            Dictionary of map values over which to calculate the processed telemetry.
            (e.g. mapDict = {'ra': [], 'dec': [], 'altitude': [], 'azimuth': [], 'airmass': []})
            Here we use 'airmass' .. which can be a single value or a numpy array.


        Returns
        -------
        numpy.ndarray, numpy.ndarray
            FWHMeff, FWHMgeom: both are the same shape numpy.ndarray.
            If airmass is a single value, FWHMeff & FWHMgeom are 1-d arrays,
            with the same order as eff_wavelen (i.e. eff_wavelen[0] = u, then FWHMeff[0] = u).
            If airmass is a numpy array, FWHMeff and FWHMgeom are 2-d arrays,
            in the order of <filter><airmass> (i.e. eff_wavelen[0] = u, 1-d array over airmass range).
        """
        fwhm_z = efdData[self.seeingcol]
        airmass = mapDict['airmass']
        airmass_correction = np.power(airmass, 0.6)
        wavelen_correction = np.power(self._config.raw_seeing_wavelength / self.eff_wavelens, 0.3)
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
