from __future__ import division
from builtins import object
from datetime import datetime
import os
import sqlite3
import numpy as np
from lsst.utils import getPackageDir

__all__ = ["SeeingData"]


class SeeingData(object):
    """Read the seeing data from disk and return appropriate FWHM_500 value at a given time.
    This is for use in simulations only. Otherwise data would come from the EFD.

    Parameters
    ----------
    time_handler : :class:`.TimeHandler`
        The instance of the simulation time handler.
    seeing_db : str or None, opt
        The name of the seeing database. If None (default), this will use the Seeing.db file
        in the 'data' directory of this package.
    """
    def __init__(self, time_handler, seeing_db=None):
        self.seeing_db = seeing_db
        if self.seeing_db is None:
            self.seeing_db = os.path.join(getPackageDir('sims_seeingModel'), 'data', 'seeing.db')
        model_time_start = datetime(time_handler.initial_dt.year, 1, 1)
        self.offset = time_handler.time_since_given_datetime(model_time_start,
                                                             reverse=True)
        self.seeing_dates = None
        self.seeing_values = None

    def __call__(self, delta_time):
        """Get the FWHM_500 value for the specified time.

        Parameters
        ----------
        delta_time : int
            The time (seconds) from the start of the simulation.

        Returns
        -------
        float
            The FWHM_500(") closest to the specified time.
        """
        delta_time += self.offset
        # Find the date to look for in the time range of the data.
        # Note that data dates should not necessarily start at zero.
        date = delta_time % self.time_range + self.min_time
        idx = np.searchsorted(self.seeing_dates, date)
        # searchsorted ensures that left < date < right
        # but we need to know if date is closer to left or to right
        left = self.seeing_dates[idx - 1]
        right = self.seeing_dates[idx]
        if date - left < right - date:
            idx -= 1
        return self.seeing_values[idx]

    def read_data(self):
        """Read the seeing information from disk.

        The default behavior is to use the module stored database. However, an
        alternate database file can be provided. The alternate database file needs to have a
        table called *Seeing* with the following columns:

        seeingId
            int : A unique index for each seeing entry.
        s_date
            int : The time (in seconds) from the start of the simulation, for the seeing observation.
        seeing
            float : The FWHM of the atmospheric PSF (in arcseconds) at zenith.
        """
        with sqlite3.connect(self.seeing_db) as conn:
            cur = conn.cursor()
            query = "select s_date, seeing from Seeing order by s_date;"
            cur.execute(query)
            results = np.array(cur.fetchall())
            self.seeing_dates = np.hsplit(results, 2)[0].flatten()
            self.seeing_values = np.hsplit(results, 2)[1].flatten()
            cur.close()
        # Make sure seeing dates are ordered appropriately (monotonically increasing).
        ordidx = self.seeing_dates.argsort()
        self.seeing_dates = self.seeing_dates[ordidx]
        self.seeing_values = self.seeing_values[ordidx]
        self.min_time = self.seeing_dates[0]
        self.max_time = self.seeing_dates[-1]
        self.time_range = self.max_time - self.min_time
