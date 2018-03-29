from .seeingData import SeeingData
from .seeingModel import SeeingModel

__all__ = ['SeeingSim']

class SeeingSim(object):
    """This is a utility class to combine SeeingData and SeeingFWHM for easier use
    with the scheduler.

    Parameters
    ----------
    time_handler : :class:`.TimeHandler`
        The instance of the simulation time handler.
    seeing_db : str or None, opt
        The name of the seeing database. If None (default), this will use the Seeing.db file
        in the 'data' directory of this package.
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
    filters_effwavelen : numpy.ndarray or None, opt
        An array containing the effective wavelengths per filter in filter_list, in nm.
        If this is None (default), sims_photUtils will be used to calculate the values for ugrizy
        based on the setup throughputs repository.
    """
    def __init__(self,  time_handler, seeing_db=None, telescope_seeing=0.25,
                 optical_design_seeing=0.08, camera_seeing=0.30,
                 raw_seeing_wavelength=500,
                 filter_effwavelens=None, filter_list=None):
        # Set up seeing data, including reading from disk.
        self.seeing_data = SeeingData(time_handler=time_handler, seeing_db=seeing_db)
        self.seeing_data.read_data()
        # Set up seeing model, including the filter wavelengths.
        self.seeing_model = SeeingModel(telescope_seeing=telescope_seeing,
                                       optical_design_seeing=optical_design_seeing,
                                       camera_seeing=camera_seeing,
                                       raw_seeing_wavelength=raw_seeing_wavelength,
                                       filter_effwavelens=filter_effwavelens)
        if filter_effwavelens is not None:
            if filter_list is None:
                raise ValueError('If filter_effwavelens is specified, so must filter_list.')
            if len(filter_list) != len(filter_effwavelens):
                raise ValueError('Length of filter_list and filter_effwavelens must match.')
            self.seeing_model.filter_list = filter_list
        self.filter_list = self.seeing_model.filter_list

    def get_fwhm500(self, delta_time):
        """Get only the FWHM500 at a given time.

        Parameters
        ----------
        delta_time : float
            The time (seconds) from the start of the simulation.

        Returns
        -------
        float
            The FWHM500.
        """
        return self.seeing_data.fwhm500_at_time(delta_time)

    def calculate_seeing(self, delta_time, filter_name, airmass):
        """Calculate seeing in a single filter -- for backwards compatibility only.

        Parameters
        ----------
        delta_time : float
            The time (seconds) from the start of the simulation.
        filter_name : str
            The single character filter name for the calculation.
        airmass : float
            The airmass for the calculation.
        Returns
        -------
        tuple
            The FWHM 500nm, FWHM Geometric and FWHM Effective seeing values.
        """
        fwhm_500 = self.seeing_data.fwhm500_at_time(delta_time)
        fwhm_eff, fwhm_geom = self.seeing_model.seeing_at_airmass(fwhm_500, airmass=airmass)
        idx = self.filter_list.index(filter_name)
        fwhm_eff = fwhm_eff[idx]
        fwhm_geom = fwhm_geom[idx]
        return (fwhm_500, fwhm_geom, fwhm_eff)

    def get_seeing(self, delta_time, airmass):
        """Calculate seeing in all filters.

        Parameters
        ----------
        delta_time : float
            The time (seconds) from the start of the simulation.
        airmass : float or np.ndarray
            The airmass value(s) at which to calculate seeing.

        Returns
        -------
        tuple
            The FWHM 500nm, FWHM Geometric and FWHM Effective seeing values,
            where FWHM_geom and FWHM_eff are np.ndarrays (1-d if airmass is float, 2-d if array).
        """
        fwhm_500 = self.seeing_data.fwhm500_at_time(delta_time)
        fwhm_eff, fwhm_geom = self.seeing_model.seeing_at_airmass(fwhm_500, airmass=airmass)
        return (fwhm_500, fwhm_geom, fwhm_eff)