import unittest
import pandas as pd
import numpy as np
from datetime import datetime
import navpy

from ephemeris_manager import EphemerisManager
from solution import func_sol


class TestGNSSProcessing(unittest.TestCase):

    def test_main(self):
        final_df = func_sol()
        self.assertFalse(final_df.empty)
        self.assertTrue("GPS_time" in final_df.columns)
        self.assertTrue("Lat" in final_df.columns)
        self.assertTrue("Alt" in final_df.columns)
        self.assertTrue("Lon" in final_df.columns)
        self.assertTrue("SatX" in final_df.columns)
        self.assertTrue("SatY" in final_df.columns)
        self.assertTrue("SatZ" in final_df.columns)


if __name__ == "__main__":
    unittest.main()
