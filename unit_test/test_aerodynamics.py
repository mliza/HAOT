import unittest
from haot.aerodynamics import *


class TestAerodynamics(unittest.TestCase):

    def setUp(self):
        """Initialize test parameters."""
        self.error_precision = 4
        self.sea_level_temperature = 300.0
        self.negative_temperature = -2.0
        self.low_mach = 0.2
        self.negative_mach = -1.0
        self.high_mach = 10.0
        self.valid_mach = 2.0

    # Test sutherland_law_viscosity #
    def test_sutherland_law_viscosity_sea_level_temperature(self):
        """Test viscosity at sea level temperature."""
        expected_viscosity = 1.568e-5
        result = sutherland_law_viscosity(self.sea_level_temperature, "Air")
        self.assertAlmostEqual(result, expected_viscosity, places=self.error_precision)

    def test_sutherland_law_viscosity_negative_temperature(self):
        """Test that a negative temperature raises a ValueError."""
        with self.assertRaises(ValueError):
            sutherland_law_viscosity(self.negative_temperature, "N2")

    def test_sutherland_law_viscosity_wrong_molecule(self):
        """Test that a wrong molecule raises a ValueError."""
        with self.assertRaises(ValueError):
            sutherland_law_viscosity(self.sea_level_temperature, "H2")

    # Test sutherland_law_viscosity #

    # Test sutherland_law_conductivity #
    def test_sutherland_law_conductivity_sea_level_temperature(self):
        """Test thermal conductivity at sea level temperature."""
        expected_viscosity = 0.0255
        result = sutherland_law_conductivity(self.sea_level_temperature, "Air")
        self.assertAlmostEqual(result, expected_viscosity, places=self.error_precision)

    def test_sutherland_law_conductivity_negative_temperature(self):
        """Test that a negative temperature raises a ValueError."""
        with self.assertRaises(ValueError):
            sutherland_law_conductivity(self.negative_temperature, "N2")

    def test_sutherland_law_conductivity_wrong_molecule(self):
        """Test that a wrong molecule raises a ValueError."""
        with self.assertRaises(ValueError):
            sutherland_law_conductivity(self.sea_level_temperature, "H2")

    # Test sutherland_law_conductivity #

    # Test speed of sound #
    def test_speed_of_sound_at_sea_level(self):
        """Test thermal conductivity at sea level temperature."""
        expected_speed_of_sound = 347.17506
        result = speed_of_sound(self.sea_level_temperature)
        self.assertAlmostEqual(
            result, expected_speed_of_sound, places=self.error_precision
        )

    def test_speed_of_sound_negative_temperature(self):
        """Test that a negative temperature raises a ValueError."""
        with self.assertRaises(ValueError):
            speed_of_sound(self.negative_temperature)

    # Test speed of sound #

    # Test isentropic relations #
    def test_isentropic_relations_standard_conditions(self):
        """Test isentropic relations for a standard Mach number."""
        mach_1 = self.valid_mach
        adiabatic_indx = 1.4
        expected_pressure_s = (1 + (adiabatic_indx - 1) / 2 * mach_1**2) ** (
            adiabatic_indx / (adiabatic_indx - 1)
        )
        expected_temperature_s = 1 + (adiabatic_indx - 1) / 2 * mach_1**2
        expected_density_s = (1 + (adiabatic_indx - 1) / 2 * mach_1**2) ** (
            1 / (adiabatic_indx - 1)
        )

        result = isentropic_relations(mach_1, adiabatic_indx)
        self.assertAlmostEqual(
            result["pressure_s"], expected_pressure_s, places=self.error_precision
        )
        self.assertAlmostEqual(
            result["temperature_s"], expected_temperature_s, places=self.error_precision
        )
        self.assertAlmostEqual(
            result["density_s"], expected_density_s, places=self.error_precision
        )

    def test_isentropic_relations_high_mach(self):
        """Test isentropic relations for a high Mach number."""
        adiabatic_indx = 1.4
        result = isentropic_relations(self.high_mach, adiabatic_indx)

        self.assertTrue(result["pressure_s"] > 1.0)
        self.assertTrue(result["temperature_s"] > 1.0)
        self.assertTrue(result["density_s"] > 1.0)

    def test_isentropic_relations_low_mach(self):
        """Test isentropic relations for a low Mach number."""
        adiabatic_indx = 1.4
        result = isentropic_relations(self.low_mach, adiabatic_indx)

        self.assertAlmostEqual(
            result["pressure_s"], 1.02828, places=self.error_precision
        )
        self.assertAlmostEqual(
            result["temperature_s"], 1.008, places=self.error_precision
        )
        self.assertAlmostEqual(result["density_s"], 1.0201, places=self.error_precision)

    def test_isentropic_relations_custom_gamma(self):
        """Test isentropic relations for a custom adiabatic index."""
        mach_1 = self.valid_mach
        adiabatic_indx = 1.3
        expected_temperature_s = 1 + (adiabatic_indx - 1) / 2 * mach_1**2
        expected_pressure_s = expected_temperature_s ** (
            adiabatic_indx / (adiabatic_indx - 1)
        )
        expected_density_s = expected_temperature_s ** (1 / (adiabatic_indx - 1))

        result = isentropic_relations(mach_1, adiabatic_indx)
        self.assertAlmostEqual(
            result["pressure_s"], expected_pressure_s, places=self.error_precision
        )
        self.assertAlmostEqual(
            result["temperature_s"], expected_temperature_s, places=self.error_precision
        )
        self.assertAlmostEqual(
            result["density_s"], expected_density_s, places=self.error_precision
        )

    def test_isentropic_relations_invalid_mach(self):
        """Test that invalid Mach numbers raise a ValueError."""
        with self.assertRaises(ValueError):
            isentropic_relations(0.0)

        with self.assertRaises(ValueError):
            isentropic_relations(self.negative_mach)

    # Test isentropic relations #


if __name__ == "__main__":
    unittest.main()
