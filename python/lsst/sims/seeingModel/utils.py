import os
import warnings


__all__ = ['get_effwavelens']

DEFAULT_THROUGHPUTS_VERSION = '1.1'
DEFAULT_FILTER_LIST = ['u', 'g', 'r', 'i', 'z', 'y']
DEFAULT_WAVELENGTHS = [367.06988658, 482.68517118,
                       622.32403587, 754.59752265,
                       869.09018708, 971.02780848]


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
    provided with the utility.
    """
    try:
        from lsst.sims.photUtils import Bandpass
        import lsst.sims.photUtils.version as photUtils_version
        no_photUtils = False
    except ImportError:
        no_photUtils = True

    fdir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    if no_photUtils or (fdir is None):
        warnings.warn('Cannot calculate effective wavelengths; either sims_photUtils is '
                      'unavailable (setup sims_photUtils) or $LSST_THROUGHPUTS_DEFAULT '
                      'is undefined (setup throughputs package). '
                      'Without these, simply using default effective wavelengths from version %s.'
                      % (DEFAULT_THROUGHPUTS_VERSION), Warning)
        photUtilsVersion = 'none'
        throughputsVersion = DEFAULT_THROUGHPUTS_VERSION
        effwavelens = []
        for f in filter_list:
            idx = filter_list.index(f)
            effwavelens.append(DEFAULT_WAVELENGTHS[idx])
    else:
        # Read the throughputs curves from the throughputs package.
        # Note that if sims_photUtils is setup, the throughputs package is as well.
        photUtilsVersion = photUtils_version.__version__
        effwavelens = []
        for f in filter_list:
            bp = Bandpass()
            bp.readThroughput(os.path.join(fdir, 'total_' + f + '.dat'))
            effwavelens.append(bp.calcEffWavelen()[1])
        with open(os.path.join(fdir, 'version_info'), 'r') as version_info:
            throughputsVersion = version_info.read().strip()
    return photUtilsVersion, throughputsVersion, effwavelens
