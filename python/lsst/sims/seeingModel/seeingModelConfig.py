import os
import warnings
import lsst.pex.config as pexConfig


__all__ = ['SeeingModelConfig', 'get_effwavelens']


DEFAULT_WAVELENGTH_VERSION = '1.3'
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


def get_effwavelens(filter_list=('u', 'g', 'r', 'i', 'z', 'y')):
    """Calculate the effective wavelengths for 'filter_list'.

    Parameters
    ----------
    filter_list: list of str, opt
        List of the filters for which to calculate effective wavelengths.
        Default = ('u', 'g', 'r', 'i', 'z', 'y')

    Returns
    -------
    List of floats
        The effective wavelengths (in nm).

    This function will attempt to calculate the effective wavelengths using
    the throughputs curves in the throughput directory and sims_photUtils.

    If sims_photUtils or throughputs is unavailable, it will just use the default values
    provided with the configuration class.
    """
    try:
        from lsst.sims.photUtils import Bandpass
        no_photUtils = False
    except ImportError:
        no_photUtils = True

    fdir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    if no_photUtils or (fdir is None):
        warnings.warn('Cannot calculate effective wavelengths; either sims_photUtils is '
                      'unavailable (setup sims_photUtils) or $LSST_THROUGHPUTS_DEFAULT '
                      'is undefined (setup throughputs package). '
                      'Without these, simply using default effective wavelengths from version %s.'
                      % (DEFAULT_WAVELENGTH_VERSION), Warning)
        effwavelens = []
        for f in filter_list:
            idx = filter_list.index(f)
            effwavelens.append(DEFAULT_WAVELENGTHS[idx])
    else:
        # Read the throughputs curves from the throughputs package.
        # Note that if sims_photUtils is setup, the throughputs package is as well.
        effwavelens = []
        for f in filter_list:
            bp = Bandpass()
            bp.readThroughput(os.path.join(fdir, 'total_' + f + '.dat'))
            effwavelens.append(bp.calcEffWavelen()[1])
    return effwavelens
