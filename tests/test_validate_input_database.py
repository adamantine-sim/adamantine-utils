# SPDX-FileCopyrightText: Copyright (c) 2026 the adamantine authors.
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from adamantine_utils.validate_input_database.validate_input_database import ValidationError, main, validate_input_database


class ValidateInputDatabaseTest(unittest.TestCase):
    def setUp(self):
        self.database = {
            "boundary": {"type": "adiabatic"},
            "discretization": {"thermal": {"fe_degree": 1}},
            "geometry": {
                "dim": 3,
                "import_mesh": False,
                "length": 1.0,
                "height": 1.0,
                "width": 1.0,
            },
            "materials": {
                "n_materials": 1,
                "property_format": "polynomial",
                "material_0": {
                    "solid": {
                        "density": 10.0,
                        "specific_heat": 10.0,
                        "thermal_conductivity_x": 10.0,
                        "thermal_conductivity_y": 10.0,
                        "thermal_conductivity_z": 10.0,
                    }
                },
            },
            "physics": {"thermal": True, "mechanical": False},
            "post_processor": {"filename_prefix": "output"},
            "refinement": {"n_heat_refinements": 0},
            "sources": {
                "n_beams": 1,
                "beam_0": {
                    "type": "goldak",
                    "scan_path_file": "scan_path.txt",
                    "scan_path_file_format": "segment",
                    "depth": 0.1,
                    "absorption_efficiency": 0.1,
                },
            },
            "time_stepping": {"method": "forward_euler", "time_step": 0.1},
        }

    def assert_invalid(self, database):
        with self.assertRaises(ValidationError):
            validate_input_database(database)

    def test_minimal_database_is_valid(self):
        validate_input_database(self.database)

    def test_optional_sections_and_conditional_properties(self):
        database = copy.deepcopy(self.database)
        database["boundary"]["type"] = "convective,traction_free"
        material = database["materials"]["material_0"]
        material["solid"]["convection_heat_transfer_coef"] = 10.0
        material["convection_temperature_infty"] = 300.0
        database["memory_space"] = "host"
        database["data_assimilation"] = {
            "convergence_tolerance": 0.01,
            "localization_cutoff_function": "gaspari_cohn",
        }
        database["units"] = {
            "mesh": "millimeter",
            "heat_source": {
                "power": "watt",
                "velocity": "meter/second",
                "dimension": "millimeter",
                "scan_path": "millimeter",
            },
        }
        database["experiment"] = {
            "read_in_experimental_data": True,
            "file": "measurements.csv",
            "format": "point_cloud",
            "first_frame": 0,
            "last_frame": 1,
            "first_camera_id": 0,
            "last_camera_id": 0,
            "log_filename": "experiment.log",
        }
        validate_input_database(database)

    def test_invalid_boundary_combination(self):
        database = copy.deepcopy(self.database)
        database["boundary"]["type"] = "adiabatic,convective"
        self.assert_invalid(database)

    def test_invalid_conditional_material_property(self):
        database = copy.deepcopy(self.database)
        database["boundary"]["type"] = "radiative"
        self.assert_invalid(database)

    def test_gaussian_beam_type_is_case_insensitive(self):
        for beam_type in ("gaussian", "GaUsSiAn"):
            with self.subTest(beam_type=beam_type):
                database = copy.deepcopy(self.database)
                database["sources"]["beam_0"]["type"] = beam_type
                validate_input_database(database)

    def test_unknown_beam_type_reports_gaussian_as_valid(self):
        database = copy.deepcopy(self.database)
        database["sources"]["beam_0"]["type"] = "unknown"
        with self.assertRaisesRegex(
            ValidationError, "Valid options are: .*'gaussian'"
        ):
            validate_input_database(database)

    def test_displacement_boundary_replaces_clamping(self):
        for dim, components in (
            (2, {"x": "t*z", "z": 0}),
            (3, {"x": "t*z", "y": "0", "z": 0}),
        ):
            with self.subTest(dim=dim):
                database = copy.deepcopy(self.database)
                database["physics"]["mechanical"] = True
                database["geometry"]["dim"] = dim
                database["boundary"] = {
                    "boundary_0": {
                        "type": "adiabatic,displacement",
                        "displacement": components,
                    },
                    "boundary_1": {"type": "adiabatic,traction_free"},
                }
                validate_input_database(database)

    def test_invalid_displacement_boundary(self):
        for section, entry, message in (
            ("type", "adiabatic,displacement", "explicit boundary"),
            ("boundary_0", {"type": "adiabatic,displacement"}, "component formula"),
            (
                "boundary_0",
                {"type": "adiabatic,displacement", "displacement": {}},
                "component formula",
            ),
            (
                "boundary_0",
                {"type": "adiabatic,displacement", "displacement": {"y": "0"}},
                "Invalid displacement component: y",
            ),
            (
                "boundary_0",
                {"type": "adiabatic,displacement", "displacement": {"x": ""}},
                "Empty displacement formula",
            ),
            (
                "boundary_0",
                {"type": "adiabatic,clamped,displacement", "displacement": {"x": "0"}},
                "cannot be combined",
            ),
            ("type", "adiabatic,clamped,clamped", "cannot be combined"),
            ("type", "adiabatic,traction_free", "needs a prescribed displacement"),
        ):
            with self.subTest(section=section, entry=entry):
                database = copy.deepcopy(self.database)
                database["physics"]["mechanical"] = True
                database["geometry"]["dim"] = 2
                database["boundary"] = {section: entry}
                with self.assertRaisesRegex(ValidationError, message):
                    validate_input_database(database)

    def test_cli_reads_json_and_returns_status(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json") as input_file:
            json.dump(self.database, input_file)
            input_file.flush()
            self.assertEqual(main([input_file.name]), 0)


if __name__ == "__main__":
    unittest.main()
