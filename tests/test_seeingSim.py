import os
import numpy as np
import unittest
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.seeingModel import SeeingSim
from lsst.sims.utils import TimeHandler


class TestSeeingSim(unittest.TestCase):
    def setUp(self):
        # Define a default set of effective wavelengths for testing.
        self.filter_list = ('u', 'g', 'r', 'i', 'z', 'y')
        self.effwavelens = np.array([ 367.06988658,  482.68517118,
                                      622.32403587,  754.59752265,
                                      869.09018708,  971.02780848])
        # Define time handler and seeing_db for data.
        self.timeHandler = TimeHandler("2020-01-01")
        # Time into seeing_db to retrieve seeing.
        self.time = 120005
        self.seeing_db = os.path.join(getPackageDir('sims_seeingModel'), 'data', 'seeing.db')
        self.seeingSim = SeeingSim(self.timeHandler)

    def test_setup(self):
        # Test various methods for setting up seeingSim.
        flist = ['u', 'g', 'r']
        feffwavelens = np.array([360.0, 485.0, 620.0])
        seeingSim = SeeingSim(self.timeHandler, seeing_db=self.seeing_db,
                              filter_list=flist,
                              filter_effwavelens=feffwavelens)
        self.assertEqual(seeingSim.filter_list, flist)
        self.assertTrue(np.all(seeingSim.seeing_model.filter_effwavelens == feffwavelens))
        self.assertEqual(seeingSim.seeing_data.seeing_db, self.seeing_db)
        self.assertEqual(len(seeingSim.seeing_data.seeing_dates), 210384)
        # Check that exception is raised appropriately if filter_list not provided when wavelens are.
        self.assertRaises(ValueError, SeeingSim, self.timeHandler, seeing_db=self.seeing_db,
                          filter_list=None, filter_effwavelens=feffwavelens)
        # Check that exception is raised appropriately if filter list is different length from wavelens.
        self.assertRaises(ValueError, SeeingSim, self.timeHandler, seeing_db=self.seeing_db,
                          filter_list=['u', 'g'], filter_effwavelens=feffwavelens)
        # Check that seeingModel sets up appropriately if not passed any information.
        seeingSim = SeeingSim(self.timeHandler)
        self.assertEqual(seeingSim.seeing_data.seeing_db, self.seeing_db)
        for i, f in enumerate(self.filter_list):
            self.assertEqual(seeingSim.filter_list[i], f)
        self.assertEqual(len(seeingSim.seeing_model.filter_effwavelens), len(seeingSim.filter_list))

    def test_get_fwhm500(self):
        # Tested full functionality in SeeingData, but here check that is passing through seeingSim ok.
        fwhm500 = self.seeingSim.get_fwhm500(self.time)
        self.assertTrue(isinstance(fwhm500, np.float))

    def test_calculate_seeing(self):
        # Check we get one filter with calculate_seeing (one airmass)
        airmass = 1.0
        fwhm500, fwhmgeom, fwhmeff = self.seeingSim.calculate_seeing(self.time, 'g', airmass)
        self.assertTrue(isinstance(fwhm500, np.float))
        self.assertTrue(isinstance(fwhmgeom, np.float))
        self.assertTrue(isinstance(fwhmeff, np.float))
        # Check that we get one filter with calculate_seeing (multiple airmasses)
        airmass = np.arange(1.0, 1.5, 0.2)
        fwhm500, fwhmgeom, fwhmeff = self.seeingSim.calculate_seeing(self.time, 'g', airmass)
        self.assertTrue(isinstance(fwhm500, np.float))
        self.assertEqual(len(fwhmgeom), len(airmass))
        self.assertEqual(len(fwhmeff), len(airmass))
        # Note that we already checked the *values* computed by these in test_seeingModel.py

    def test_get_seeing(self):
        # Check that we get multiple filters with get_seeing (one airmass)
        airmass = 1.0
        fwhm500, fwhmgeom, fwhmeff = self.seeingSim.get_seeing(self.time, airmass)
        self.assertTrue(isinstance(fwhm500, np.float))
        self.assertEqual(fwhmgeom.shape, (len(self.filter_list), ))
        self.assertEqual(fwhmeff.shape, (len(self.filter_list), ))
        # Check that we get one filter with calculate_seeing (multiple airmasses)
        airmass = np.arange(1.0, 1.5, 0.2)
        fwhm500, fwhmgeom, fwhmeff = self.seeingSim.get_seeing(self.time, airmass)
        self.assertTrue(isinstance(fwhm500, np.float))
        self.assertEqual(fwhmgeom.shape, (len(self.filter_list), len(airmass)))
        self.assertEqual(fwhmeff.shape, (len(self.filter_list), len(airmass)))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()