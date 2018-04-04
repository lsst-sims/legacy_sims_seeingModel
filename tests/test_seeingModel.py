import numpy as np
import unittest
import lsst.utils.tests
from lsst.sims.seeingModel import SeeingModel


class TestSeeingModel(unittest.TestCase):
    def setUp(self):
        # Define a default set of effective wavelengths for testing.
        self.effwavelens = np.array([ 367.06988658,  482.68517118,
                                      622.32403587,  754.59752265,
                                      869.09018708,  971.02780848])

    def test_fwhm_system_zenith(self):
        # Check calculation is being done as expected.
        seeingModel = SeeingModel(telescope_seeing=0.25,
                                  optical_design_seeing=0.08, camera_seeing=0.30,
                                  filter_effwavelens=self.effwavelens)
        self.assertEqual(seeingModel.fwhm_system_zenith, 0.39862262855989494)
        # Check defaults haven't changed unexpectedly.
        seeingModel = SeeingModel(filter_effwavelens=self.effwavelens)
        self.assertEqual(seeingModel.fwhm_system_zenith, 0.39862262855989494)

    def test_set_effwavelens(self):
        # Check the values are set as expected when we set them.
        seeingModel = SeeingModel(filter_effwavelens=self.effwavelens)
        self.assertTrue(np.all(self.effwavelens==seeingModel.filter_effwavelens))
        # Check that the values can be set using sims_photUtils.
        seeingModel = SeeingModel()
        self.assertTrue(seeingModel.filter_effwavelens is not None)

    def test_fwhmGeomEff(self):
        # Check that the translation between FWHM effective and geometric is done as expected.
        # (note that fwhmEff_tofwhmGeom & fwhmGeom_to_fwhmEff are static methods)
        # Document-20160 for reference.
        fwhm_eff = 1.23
        fwhm_geom = 0.822 * fwhm_eff + 0.052
        self.assertEqual(fwhm_geom, SeeingModel.fwhmEff_to_fwhmGeom(fwhm_eff))
        self.assertEqual(fwhm_eff, SeeingModel.fwhmGeom_to_fwhmEff(fwhm_geom))

    def test_seeingCalc(self):
        # Check the calculation from fwhm_500 to fwhm_eff/fwhm_geom.
        # Use simple effective wavelengths and airmass values.
        wavelens = np.array([500.0, 1000.0])
        seeingModel = SeeingModel(telescope_seeing=0.0,
                                  optical_design_seeing=0.0, camera_seeing=0.0,
                                  filter_effwavelens=wavelens)
        fwhm_500 = 1.0
        # Single airmass.
        airmass = 1.0
        fwhm_eff, fwhm_geom = seeingModel.seeing_at_airmass(fwhm_500, airmass=airmass)
        self.assertEqual(len(fwhm_geom), len(seeingModel.filter_effwavelens))
        self.assertEqual(fwhm_eff.shape, (len(seeingModel.filter_effwavelens),))
        # Check seeing in u (see filter_list in seeingModel) for order
        expected_fwhm_eff = fwhm_500 * 1.16 * np.sqrt(1.04)
        self.assertEqual(fwhm_eff[0], expected_fwhm_eff)
        # Check scaling with wavelength.
        expected_fwhm_eff = 1.16 * np.sqrt(1.04) * fwhm_500 * np.power(500./wavelens[1], 0.3)
        self.assertAlmostEqual(fwhm_eff[1], expected_fwhm_eff, places=15)
        # Multiple airmasses.
        airmass = np.array([1.0, 1.5])
        fwhm_eff, fwhm_geom = seeingModel.seeing_at_airmass(fwhm_500, airmass=airmass)
        self.assertEqual(fwhm_eff.shape, (len(seeingModel.filter_effwavelens),len(airmass)))
        expected_fwhm_eff = fwhm_500 * 1.16 * np.sqrt(1.04)
        self.assertEqual(fwhm_eff[0][0], expected_fwhm_eff)
        # Check scaling with airmass.
        expected_fwhm_eff = expected_fwhm_eff * np.power(airmass[1], 0.6)
        self.assertAlmostEqual(fwhm_eff[0][1], expected_fwhm_eff, places=15)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
