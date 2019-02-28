import numpy as np
import unittest
import lsst.utils.tests
from lsst.sims.seeingModel import SeeingModel, SeeingModelConfig


class TestSeeingModel(unittest.TestCase):
    def setUp(self):
        # Set config to known values.
        config = SeeingModelConfig()
        config.telescope_seeing = 0.25
        config.optical_design_seeing = 0.08
        config.camera_seeing = 0.30
        config.raw_seeing_wavelength = 500
        # Define a non-default set of effective wavelengths for testing.
        filterlist = ['u', 'g', 'r', 'i', 'z']
        effwavelens = [ 367.06988658, 482.68517118,
                        622.32403587,  754.59752265,
                        869.09018708]
        config.filter_list = filterlist
        config.filter_effwavelens = effwavelens
        config.efd_columns = ['raw_seeing']
        config.efd_delta_time = 0
        self.config = config

    def test_configure(self):
        # Configure with defaults.
        seeingModel = SeeingModel()
        seeingModel.configure()
        conf = SeeingModelConfig()
        self.assertEqual(seeingModel.config, conf)
        # Test specifying the config.
        seeingModel = SeeingModel()
        seeingModel.configure(self.config)
        self.assertEqual(seeingModel.config.filter_list, self.config.filter_list)
        # Test specifying an incorrect config.
        self.assertRaises(ValueError, seeingModel.configure, 0.8)
        # Test specifying a config override filename.
        # TBD.

    def test_fwhm_system_zenith(self):
        # Check calculation is being done as expected.
        seeingModel = SeeingModel()
        seeingModel.configure()
        self.assertEqual(seeingModel.fwhm_system_zenith, 0.39862262855989494)

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
        seeingModel = SeeingModel()
        seeingModel.configure(self.config)
        # Hack the effective wavelengths in the class.
        seeingModel.eff_wavelens = wavelens
        # Simple fwhm_500 input.
        fwhm_500 = 1.0
        # Single airmass.
        airmass = 1.0
        fwhm_eff, fwhm_geom = seeingModel.seeing_at_airmass(fwhm_500, airmass=airmass)
        # Check shape of returned values.
        self.assertEqual(fwhm_eff.shape, (len(seeingModel.eff_wavelens),))
        # Check actual value of seeing in @ wavelen[0] @ zenith after addition of system.
        fwhm_system = seeingModel.fwhm_system_zenith
        expected_fwhm_eff = 1.16 * np.sqrt(fwhm_system ** 2 + 1.04 * fwhm_500 ** 2)
        self.assertAlmostEqual(fwhm_eff[0], expected_fwhm_eff, 15)
        # Check expected value if we remove the system component.
        seeingModel.fwhm_system_zenith = 0
        fwhm_eff, fwhm_geom = seeingModel.seeing_at_airmass(fwhm_500, airmass=airmass)
        expected_fwhm_eff = 1.16 * np.sqrt(1.04) * fwhm_500
        self.assertAlmostEqual(fwhm_eff[0], expected_fwhm_eff, 15)
        # Check scaling with wavelength (remove system component).
        expected_fwhm_eff = 1.16 * np.sqrt(1.04) * fwhm_500 * np.power(500./wavelens[1], 0.3)
        self.assertAlmostEqual(fwhm_eff[1], expected_fwhm_eff, places=15)
        # Multiple airmasses.
        airmass = np.array([1.0, 1.5])
        fwhm_eff, fwhm_geom = seeingModel.seeing_at_airmass(fwhm_500, airmass=airmass)
        self.assertEqual(fwhm_eff.shape, (len(seeingModel.eff_wavelens), len(airmass)))
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
