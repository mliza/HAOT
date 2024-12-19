import unittest
from haot.quantum_mechanics import *


class TestOptics(unittest.TestCase):

    def setUp(self):
        """Initialize test parameters."""
        self.error_precision = 4
        self.valid_temperature = 200.0
        self.invalid_temperature = -20.0
        self.valid_molecule = "H2"
        self.invalid_molecule = "Argon"
        self.valid_wavelength = 633.0
        self.invalid_wavelength = -2300.0

    # Test zero_point_energy #
    def test_zero_point_energy_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            zero_point_energy(self.invalid_molecule)
    # Test zero_point_energy #

if __name__ == "__main__":
    unittest.main()
