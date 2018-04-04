import os
import unittest
import sqlite3
from lsst.utils import getPackageDir
import lsst.utils.tests
from lsst.utils.tests import getTempFilePath
from lsst.sims.seeingModel import SeeingData
from lsst.sims.utils import TimeHandler


class TestSeeingData(unittest.TestCase):

    def setUp(self):
        self.timeHandler = TimeHandler("2020-01-01")
        self.seeing_db = os.path.join(getPackageDir('sims_seeingModel'), 'data', 'seeing.db')

    def test_basic_information_after_creation(self):
        seeingData = SeeingData(self.timeHandler, seeing_db=self.seeing_db)
        self.assertIsNone(seeingData.seeing_dates)
        self.assertIsNone(seeingData.seeing_values)
        self.assertEqual(seeingData.offset, 0)
        self.assertEqual(seeingData.seeing_db, self.seeing_db)
        # And check sets seeing_db appropriately if not provided.
        seeingData = SeeingData(self.timeHandler, seeing_db=None)
        self.assertEqual(seeingData.seeing_db, self.seeing_db)

    def test_information_after_read(self):
        seeingData = SeeingData(self.timeHandler, seeing_db=self.seeing_db)
        seeingData.read_data()
        self.assertTrue(len(seeingData.seeing_values) > 0)
        self.assertTrue(len(seeingData.seeing_dates) > 0)

    def test_fwhm500_at_time(self):
        seeingData = SeeingData(self.timeHandler, self.seeing_db)
        seeingData.read_data()
        self.assertEqual(seeingData.fwhm500_at_time(75400), 0.859431982040405)
        self.assertEqual(seeingData.fwhm500_at_time(76700), 0.646009027957916)
        self.assertEqual(seeingData.fwhm500_at_time(63190400), 0.64860999584198)
        self.assertEqual(seeingData.fwhm500_at_time(189424900), 0.699440002441406)
        # Test time selection from seeing data.
        fwhm500 =  seeingData.fwhm500_at_time(800)
        # Hack seeing data to remove first date, thus db does not start at zero.
        seeingData.seeing_dates = seeingData.seeing_dates[:-1]
        seeingData.seeing_values = seeingData.seeing_values[:-1]
        seeingData.time_range = seeingData.seeing_dates[-1] - seeingData.seeing_dates[0]
        seeingData.min_time = seeingData.seeing_dates[0]
        self.assertEqual(fwhm500, seeingData.fwhm500_at_time(800))

    def test_using_different_start_month(self):
        seeingData = SeeingData(TimeHandler("2020-05-24"), self.seeing_db)
        self.assertEqual(seeingData.offset, 12441600)
        seeingData.read_data()
        self.assertEqual(seeingData.fwhm500_at_time(75400), 0.437314003705978)
        self.assertEqual(seeingData.fwhm500_at_time(76700), 0.510206997394562)
        self.assertEqual(seeingData.fwhm500_at_time(63190400), 0.453994989395142)
        self.assertEqual(seeingData.fwhm500_at_time(189424900), 0.386815994977951)

    def test_alternate_db(self):
        # Create temporary data file, use it as the seeing_db.
        with getTempFilePath('.alt_seeing.db') as tmpdb:
            seeing_table = []
            seeing_table.append("seeingId INTEGER PRIMARY KEY")
            seeing_table.append("s_date INTEGER")
            seeing_table.append("seeing DOUBLE")
            with sqlite3.connect(tmpdb) as conn:
                cur = conn.cursor()
                cur.execute("DROP TABLE IF EXISTS Seeing")
                cur.execute("CREATE TABLE Seeing({})".format(",".join(seeing_table)))
                cur.executemany("INSERT INTO Seeing VALUES(?, ?, ?)", [(1, 9997, 0.5), (2, 10342, 0.3)])
                cur.close()
            seeingData = SeeingData(self.timeHandler, seeing_db=tmpdb)
            seeingData.read_data()
            self.assertEqual(seeingData.seeing_values.size, 2)
            self.assertEqual(seeingData.seeing_values[1], 0.3)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
