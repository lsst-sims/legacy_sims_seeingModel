import lsst.pex.config as pexConfig


__all__ = ['SeeingModelConfig']

DEFAULT_THROUGHPUTS_VERSION = '1.1'
DEFAULT_FILTER_LIST = ['u', 'g', 'r', 'i', 'z', 'y']
DEFAULT_WAVELENGTHS = [367.06988658, 482.68517118,
                       622.32403587, 754.59752265,
                       869.09018708, 971.02780848]


class SeeingModelConfig(pexConfig.Config):
    """A pex_config configuration class for default seeing model parameters.
    """
    telescope_seeing = pexConfig.Field(doc="Telescope contribution to IQ (arcsec)",
                                       dtype=float,
                                       default=0.25)
    optical_design_seeing = pexConfig.Field(doc="Optics contribution to IQ (arcsec)",
                                            dtype=float,
                                            default=0.08)
    camera_seeing = pexConfig.Field(doc="Camera contribution to IQ (arcsec)",
                                    dtype=float,
                                    default=0.30)
    raw_seeing_wavelength = pexConfig.Field(doc="Wavelength of input zenith IQ (nm)",
                                            dtype=float,
                                            default=500)

    filter_list = pexConfig.ListField(doc="List of filters for which to calculate seeing",
                                      dtype=str,
                                      default=DEFAULT_FILTER_LIST)
    filter_effwavelens = pexConfig.ListField(doc="Effective wavelengths for filters (nm)",
                                             dtype=float,
                                             default=DEFAULT_WAVELENGTHS)
    throughputs_version = pexConfig.Field(doc="Version of the throughputs files",
                                          dtype=str,
                                          default=DEFAULT_THROUGHPUTS_VERSION)

    efd_columns = pexConfig.ListField(doc="List of data required from EFD",
                                      dtype=str,
                                      default=['FWHM_500'])
    efd_delta_time = pexConfig.Field(doc="Length (delta time) of history to request from the EFD (seconds)",
                                     dtype=float,
                                     default=0)
    target_columns = pexConfig.ListField(doc="Names of the keys for the airmass in the "
                                             "scheduler target maps",
                                         dtype=str,
                                         default=['airmass'])
