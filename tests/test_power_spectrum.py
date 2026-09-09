import unittest

import lightkurve as lk
import numpy as np

from asteroseismologie import compute_power_spectrum


class PowerSpectrumTests(unittest.TestCase):
    def test_parseval_normalization(self):
        count = 4096
        cadence_s = 60.0
        time_days = np.arange(count) * cadence_s / 86400.0
        flux_ppm = 100.0 * np.sin(2.0 * np.pi * 1000e-6 * np.arange(count) * cadence_s)
        light_curve = lk.LightCurve(time=time_days, flux=flux_ppm)

        freq, power = compute_power_spectrum(
            light_curve,
            fmin_uHz=1.0,
            fmax_uHz=2000.0,
            oversample=1,
            auto_raise_sc_fmax=False,
        )

        integrated_power = np.trapezoid(power, freq)
        self.assertAlmostEqual(integrated_power, np.var(flux_ppm), delta=0.01 * np.var(flux_ppm))


if __name__ == "__main__":
    unittest.main()