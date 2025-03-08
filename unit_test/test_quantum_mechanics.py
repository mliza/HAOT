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
        self.valid_vibrational = 2
        self.invalid_vibrational = -2

    def test_zero_point_energy_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            zero_point_energy(self.invalid_molecule)

    def test_vibrational_partition_function_invalid_vibrational(self):
        """Test invalid vibrational quantum number."""
        with self.assertRaises(ValueError):
            vibrational_partition_function(
                self.invalid_vibrational, self.valid_temperature, self.valid_molecule
            )

    def test_vibrational_partition_function_invalid_temperature(self):
        """Test invalid temperature."""
        with self.assertRaises(ValueError):
            vibrational_partition_function(
                self.valid_vibrational, self.invalid_temperature, self.valid_molecule
            )

    def test_vibrational_partition_function_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            vibrational_partition_function(
                self.valid_vibrational, self.valid_temperature, self.invalid_molecule
            )

    def test_rotational_partition_function_invalid_vibrational(self):
        """Test invalid rotational quantum number."""
        with self.assertRaises(ValueError):
            rotational_partition_function(
                self.invalid_rotational, self.valid_temperature, self.valid_molecule
            )

    def test_rotational_partition_function_invalid_temperature(self):
        """Test invalid temperature."""
        with self.assertRaises(ValueError):
            rotational_partition_function(
                self.valid_rotational, self.invalid_temperature, self.valid_molecule
            )

    def test_rotational_partition_function_invalid_molecule(self):
        """Test invalid molecule."""
        with self.assertRaises(ValueError):
            rotational_partition_function(
                self.valid_rotational, self.valid_temperature, self.invalid_molecule
            )


if __name__ == "__main__":
    unittest.main()
