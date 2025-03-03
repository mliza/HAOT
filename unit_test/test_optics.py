import unittest
from haot.optics import *


class TestOptics(unittest.TestCase):

    def setUp(self):
        """Initialize test parameters."""
        self.error_precision = 4
        self.valid_temperature = 200.0
        self.mass_density = 0.5
        self.invalid_temperature = -20.0
        self.valid_molecule = "H2"
        self.invalid_molecule = "Argon"
        self.valid_wavelength = 633.0
        self.invalid_wavelength = -2300.0

    # Test kerl_polarizability_temperature #
    def test_kerl_polarizability_valid_temperature(self):
        """Test invalid temperature."""
        result = kerl_polarizability_temperature(
            self.valid_temperature, self.valid_molecule, self.valid_wavelength
        )

    def test_kerl_polarizability_invalid_temperature(self):
        """Test invalid temperature."""
        with self.assertRaises(ValueError):
            kerl_polarizability_temperature(
                self.invalid_temperature, self.valid_molecule, self.valid_wavelength
            )

    def test_kerl_polarizability_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            kerl_polarizability_temperature(
                self.valid_temperature, self.invalid_molecule, self.valid_wavelength
            )

    def test_kerl_polarizability_invalid_wavelength(self):
        """Test invalid wavelength."""
        with self.assertRaises(ValueError):
            kerl_polarizability_temperature(
                self.valid_temperature, self.valid_molecule, self.invalid_wavelength
            )

    def test_index_of_refraction_density_temperature_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            index_of_refraction_density_temperature(self.valid_temperature,
                                                    self.mass_density,
                                                    self.invalid_molecule,
                                                    self.valid_wavelength
            )


if __name__ == "__main__":
    unittest.main()
