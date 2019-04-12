import unittest
import lsst.utils.tests
from lsst.sims.seeingModel import SeeingModelConfig
from lsst.sims.seeingModel import get_effwavelens


class TestSeeingModel(unittest.TestCase):
    def testConfig(self):
        # Just have to test that there are no errors in normal function.
        config = SeeingModelConfig()
        config.validate()
        config.freeze()


class TestGetEffWavelens(unittest.TestCase):
    def test_get_effwavelens(self):
        # In default setup, should get all filters.
        filters = ('u', 'g', 'r', 'i', 'z', 'y')
        photUtils_version, throughputs_version, effwavelens = get_effwavelens()
        self.assertEqual(len(filters), len(effwavelens))
        filters = ('u', 'g')
        photUtils_version, throughputs_version, effwavelens = get_effwavelens(filters)
        self.assertEqual(len(filters), len(effwavelens))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
