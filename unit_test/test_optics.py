import unittest
from haot.optics import *


class TestOptics(unittest.TestCase):

    def setUp(self):
        """Initialize test parameters."""
        self.error_precision = 4
        self.valid_temperature = 200.0
        self.mass_density = 0.5
        self.invalid_temperature = -20.0
        self.invalid_index = -2.5
        self.valid_index = np.array([2.5, 2.0])
        self.valid_molecule = "H2"
        self.invalid_molecule = "Argon"
        self.valid_wavelength = 633.0
        self.invalid_wavelength = -2300.0
        self.invalid_mass_density_dict = {"N2": 0.5, "H2": 0.1}
        self.invalid_distance = -2
        self.invalid_distance_dimension = np.array([2, 1, -4])
        self.valid_distance = 2.0

    # Test permittivity_material, invalid data type
    def test_permittivity_material_invalid_data_type(self):
        """Test invalid index of refraction data type."""
        with self.assertRaises(ValueError):
            permittivity_material(self.valid_molecule)

    # Test permittivity_material, invalid index
    def test_permittivity_material_invalid_index(self):
        """Test invalid index of refraction."""
        with self.assertRaises(ValueError):
            permittivity_material(self.invalid_index)

    # Test electric_susceptibility, invalid data type
    def test_electric_susceptibility_invalid_data_type(self):
        """Test invalid index of refraction data type."""
        with self.assertRaises(ValueError):
            electric_susceptibility(self.valid_molecule)

    # Test electric_susceptibility, invalid index
    def test_electric_susceptibility_invalid_index(self):
        """Test invalid index of refraction."""
        with self.assertRaises(ValueError):
            electric_susceptibility(self.invalid_index)

    # Test optical_path_length, invalid data type
    def test_optical_path_length_invalid_data_type(self):
        """Test invalid index of refraction data type."""
        with self.assertRaises(ValueError):
            optical_path_length(self.valid_molecule, self.valid_distance)

    # Test optical_path_length, invalid index
    def test_optical_path_length_invalid_index(self):
        """Test invalid index of refraction."""
        with self.assertRaises(ValueError):
            optical_path_length(self.invalid_index, self.valid_distance)

    # Test optical_path_length, invalid index
    def test_optical_path_length_invalid_distance(self):
        """Test invalid distance."""
        with self.assertRaises(ValueError):
            optical_path_length(self.invalid_index, self.invalid_distance)

    # Test optical_path_length, invalid index
    def test_optical_path_length_invalid_index(self):
        """Test invalid index of refraction."""
        with self.assertRaises(ValueError):
            optical_path_length(self.invalid_index, self.valid_distance)

    # Test optical_path_length, invalid index
    def test_optical_path_length_invalid_distance(self):
        """Test invalid distance."""
        with self.assertRaises(ValueError):
            optical_path_length(self.valid_index, self.invalid_distance)

    # Test optical_path_length, invalid index
    def test_optical_path_length_invalid_distance_dimension(self):
        """Test invalid dimensions."""
        with self.assertRaises(ValueError):
            optical_path_length(self.valid_index, self.invalid_distance_dimension)

    # Test index_of_refraction invalid mass density format
    def test_kerl_index_of_refraction_invalid_input(self):
        """Test invalid mass density."""
        with self.assertRaises(ValueError):
            index_of_refraction(self.mass_density)

    # Test index_of_refraction wrong keys
    def test_kerl_index_of_refraction_invalid_keys(self):
        """Test invalid mass mass density keys."""
        with self.assertRaises(ValueError):
            index_of_refraction(self.invalid_mass_density_dict)

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
            index_of_refraction_density_temperature(
                self.valid_temperature,
                self.mass_density,
                self.invalid_molecule,
                self.valid_wavelength,
            )

    def test_index_of_refraction_density_temperature(self):
        """Test invalid temperature."""
        with self.assertRaises(ValueError):
            index_of_refraction_density_temperature(
                self.invalid_temperature,
                self.mass_density,
                self.valid_molecule,
                self.valid_wavelength,
            )

    def test_index_of_refraction_density_temperature(self):
        """Test invalid wavelength."""
        with self.assertRaises(ValueError):
            index_of_refraction_density_temperature(
                self.valid_temperature,
                self.mass_density,
                self.valid_molecule,
                self.invalid_wavelength,
            )


if __name__ == "__main__":
    unittest.main()
