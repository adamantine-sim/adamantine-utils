"""Validate an Adamantine JSON input file without building Adamantine.

Usage:
    python3 validate_input_database.py path/to/input.json

The checks in this module intentionally follow
``source/validate_input_database.cc``.  Material-property table and polynomial
contents are not checked here; those checks remain the responsibility of the
corresponding Adamantine material-property implementation.
"""

# SPDX-FileCopyrightText: Copyright (c) 2026 the adamantine authors.
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional


class ValidationError(ValueError):
    """An input-file validation error."""


_MISSING = object()
_OPTIONAL_MISSING = object()


def _get(data: Any, path: str, default: Any = _MISSING) -> Any:
    """Return a dotted-path value, or raise a useful missing-value error."""

    current = data
    for component in path.split("."):
        if not isinstance(current, dict) or component not in current:
            if default is not _MISSING:
                return default
            raise ValidationError(f"Missing required input '{path}'.")
        current = current[component]
    return current


def _require_object(data: Any, path: str) -> Dict[str, Any]:
    value = _get(data, path)
    if not isinstance(value, dict):
        raise ValidationError(f"Input '{path}' must be a JSON object.")
    return value


def _as_bool(value: Any, path: str) -> bool:
    if isinstance(value, bool):
        return value
    # Accept the textual forms produced by property-tree inputs as well as
    # native JSON booleans.  This is useful when a JSON file was generated from
    # an INFO file and all scalar values are strings.
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    raise ValidationError(f"Input '{path}' must be a boolean.")


def _as_unsigned(value: Any, path: str) -> int:
    if isinstance(value, bool):
        raise ValidationError(f"Input '{path}' must be a non-negative integer.")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str):
        try:
            result = int(value.strip(), 10)
        except ValueError as error:
            raise ValidationError(
                f"Input '{path}' must be a non-negative integer."
            ) from error
    else:
        raise ValidationError(f"Input '{path}' must be a non-negative integer.")

    # unsigned int is the type used by the C++ validator.
    if result < 0 or result > 2**32 - 1:
        raise ValidationError(f"Input '{path}' must be a non-negative integer.")
    return result


def _as_number(value: Any, path: str) -> float:
    if isinstance(value, bool):
        raise ValidationError(f"Input '{path}' must be a number.")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"Input '{path}' must be a number.") from error
    if not math.isfinite(result):
        raise ValidationError(f"Input '{path}' must be a finite number.")
    return result


def _as_string(value: Any, path: str) -> str:
    if not isinstance(value, str):
        raise ValidationError(f"Input '{path}' must be a string.")
    return value


def _require_keys(data: Dict[str, Any], keys: Iterable[str], message: str) -> None:
    for key in keys:
        if key not in data:
            raise ValidationError(message)


def _validate_boundary(
    database: Dict[str, Any], use_thermal_physics: bool, use_mechanical_physics: bool
) -> tuple[bool, bool]:
    """Validate boundaries and return (uses_convective, uses_radiative)."""

    boundary_database = _require_object(database, "boundary")
    uses_convective = False
    uses_radiative = False
    mechanical_dirichlet_bc = False

    for boundary_name, boundary_entry in boundary_database.items():
        line_types = set()
        if boundary_name == "type":
            boundary_type_value = boundary_entry
            boundary_path = "boundary.type"
        else:
            boundary_entry = (
                boundary_entry
                if isinstance(boundary_entry, dict)
                else {}
            )
            if "type" not in boundary_entry:
                raise ValidationError(
                    f"Missing required input 'boundary.{boundary_name}.type'."
                )
            boundary_type_value = boundary_entry["type"]
            boundary_path = f"boundary.{boundary_name}.type"

        boundary_type_string = _as_string(boundary_type_value, boundary_path)
        thermal_bc = False
        mechanical_bc = False

        for boundary_type in boundary_type_string.split(","):
            if boundary_type == "adiabatic":
                if "radiative" in line_types and "convective" in line_types:
                    raise ValidationError(
                        "adiabatic condition cannot be combined with another "
                        "type of thermal boundary condition."
                    )
                line_types.add("adiabatic")
                thermal_bc = True
            elif boundary_type == "radiative":
                if "adiabatic" in line_types:
                    raise ValidationError(
                        "adiabatic condition cannot be combined with another "
                        "type of thermal boundary condition."
                    )
                line_types.add("radiative")
                uses_radiative = True
                thermal_bc = True
            elif boundary_type == "convective":
                if "adiabatic" in line_types:
                    raise ValidationError(
                        "adiabatic condition cannot be combined with another "
                        "type of thermal boundary condition."
                    )
                line_types.add("convective")
                uses_convective = True
                thermal_bc = True
            elif boundary_type in {"clamped", "traction_free", "displacement"}:
                if mechanical_bc:
                    raise ValidationError(
                        "Mechanical boundary conditions cannot be combined "
                        "together."
                    )
                line_types.add(boundary_type)
                mechanical_bc = True
            else:
                raise ValidationError("Unknown boundary type.")

        if "displacement" in line_types:
            if not boundary_name.startswith("boundary_"):
                raise ValidationError(
                    "Displacement requires an explicit boundary_<id> section."
                )
            components = boundary_entry.get("displacement")
            if not isinstance(components, dict) or not components:
                raise ValidationError("Displacement requires a component formula.")
            dimension = _as_unsigned(_get(database, "geometry.dim"), "geometry.dim")
            valid_components = {"x", "y", "z"} if dimension == 3 else {"x", "z"}
            for component, formula in components.items():
                if component not in valid_components:
                    raise ValidationError(f"Invalid displacement component: {component}")
                # Boost Property Tree treats JSON scalars as leaf text.
                if isinstance(formula, (dict, list)) or formula == "":
                    raise ValidationError("Empty displacement formula.")

        if line_types & {"clamped", "displacement"}:
            mechanical_dirichlet_bc = True
        if use_thermal_physics and not thermal_bc:
            raise ValidationError("Missing thermal boundary condition.")
        if use_mechanical_physics and not mechanical_bc:
            raise ValidationError("Missing mechanical boundary condition.")

    if use_mechanical_physics and not mechanical_dirichlet_bc:
        raise ValidationError(
            "At least one boundary needs a prescribed displacement or clamping."
        )

    return uses_convective, uses_radiative


def validate_input_database(database: Dict[str, Any]) -> None:
    """Validate an Adamantine JSON database.

    The function returns ``None`` when the input is valid and raises
    :class:`ValidationError` on the first invalid condition, matching the
    fail-fast behavior of the C++ validator.
    """

    if not isinstance(database, dict):
        raise ValidationError("The top-level JSON value must be an object.")

    # Tree: physics
    use_thermal_physics = _as_bool(_get(database, "physics.thermal"), "physics.thermal")
    use_mechanical_physics = _as_bool(
        _get(database, "physics.mechanical"), "physics.mechanical"
    )
    if not (use_thermal_physics or use_mechanical_physics):
        raise ValidationError("Both thermal and mechanical physics are disabled")

    # Tree: boundary
    uses_convective, uses_radiative = _validate_boundary(
        database, use_thermal_physics, use_mechanical_physics
    )

    # Tree: discretization.thermal
    if use_thermal_physics:
        fe_degree = _as_unsigned(
            _get(database, "discretization.thermal.fe_degree"),
            "discretization.thermal.fe_degree",
        )
        if not (0 < fe_degree < 6):
            raise ValidationError("fe_degree should be between 1 and 5.")

        quadrature_type = _get(
            database, "discretization.thermal.quadrature", _OPTIONAL_MISSING
        )
        if quadrature_type is not _OPTIONAL_MISSING:
            quadrature_type = _as_string(
                quadrature_type, "discretization.thermal.quadrature"
            )
            if quadrature_type.lower() not in {"gauss", "lobatto"}:
                raise ValidationError("Unknown quadrature type.")

    # Tree: geometry
    dim = _as_unsigned(_get(database, "geometry.dim"), "geometry.dim")
    if dim not in {2, 3}:
        raise ValidationError("dim should be 2 or 3")

    use_powder = _as_bool(
        _get(database, "geometry.use_powder", False), "geometry.use_powder"
    )
    if use_powder:
        powder_layer = _as_number(
            _get(database, "geometry.powder_layer"), "geometry.powder_layer"
        )
        if powder_layer < 0.0:
            raise ValidationError("powder_layer must be non-negative.")

    material_deposition = _as_bool(
        _get(database, "geometry.material_deposition", False),
        "geometry.material_deposition",
    )
    if material_deposition:
        method = _as_string(
            _get(database, "geometry.material_deposition_method"),
            "geometry.material_deposition_method",
        )
        if method not in {"file", "scan_paths"}:
            raise ValidationError(
                "Method type for material deposition, "
                f"'{method}', is not recognized. Valid options are: 'file' "
                "and 'scan_paths'"
            )
        if method == "file":
            if "material_deposition_file" not in _require_object(
                database, "geometry"
            ):
                raise ValidationError(
                    "If the material deposition method is 'file', "
                    "'material_deposition_file' must be given."
                )
        else:
            geometry = _require_object(database, "geometry")
            _require_keys(
                geometry,
                ["deposition_length"],
                "If the material deposition method is 'scan_path', "
                "'deposition_length' must be given.",
            )
            _require_keys(
                geometry,
                ["deposition_height"],
                "If the material deposition method is 'scan_path', "
                "'deposition_height' must be given.",
            )
            if dim == 3:
                _require_keys(
                    geometry,
                    ["deposition_width"],
                    "If the material deposition method is 'scan_path', "
                    "'deposition_width' must be given.",
                )
            _require_keys(
                geometry,
                ["deposition_lead_time"],
                "If the material deposition method is 'scan_path', "
                "'deposition_lead_time' must be given.",
            )

    import_mesh = _as_bool(
        _get(database, "geometry.import_mesh"), "geometry.import_mesh"
    )
    geometry = _require_object(database, "geometry")
    if import_mesh:
        _require_keys(
            geometry,
            ["mesh_file"],
            "If the the mesh is imported, 'mesh_file' must be given.",
        )
        _require_keys(
            geometry,
            ["mesh_format"],
            "If the the mesh is imported, 'mesh_format' must be given.",
        )
    else:
        _require_keys(
            geometry,
            ["length"],
            "If the the mesh is not imported, 'length' must be given.",
        )
        _require_keys(
            geometry,
            ["height"],
            "If the the mesh is not imported, 'height' must be given.",
        )
        if dim == 3:
            _require_keys(
                geometry,
                ["width"],
                "If the the mesh is not imported, 'width' must be given.",
            )

    # Tree: materials
    materials = _require_object(database, "materials")
    n_materials = _as_unsigned(
        _get(database, "materials.n_materials"), "materials.n_materials"
    )
    property_format = _as_string(
        _get(database, "materials.property_format"),
        "materials.property_format",
    )
    if property_format not in {"table", "polynomial"}:
        raise ValidationError("property_format should be table or polynomial.")

    for material_id in range(n_materials):
        material_name = f"material_{material_id}"
        if material_name not in materials:
            raise ValidationError(
                "Number of material subtrees does not match the set number of "
                "materials."
            )
        material = materials[material_name]
        if not isinstance(material, dict):
            raise ValidationError(f"Input 'materials.{material_name}' must be a JSON object.")

        has_a_valid_state = False
        for state_name in ("solid", "liquid", "powder"):
            if state_name not in material:
                continue
            has_a_valid_state = True
            state = material[state_name]
            if not isinstance(state, dict):
                raise ValidationError(
                    f"Input 'materials.{material_name}.{state_name}' must be a JSON object."
                )

            _require_keys(
                state,
                ["density"],
                "Each state needs a user-specified density.",
            )
            _require_keys(
                state,
                ["specific_heat"],
                "Each state needs a user-specified specific heat.",
            )
            _require_keys(
                state,
                ["thermal_conductivity_x"],
                "Each state needs a user-specified specific thermal conductivity x.",
            )
            _require_keys(
                state,
                ["thermal_conductivity_z"],
                "Each state needs a user-specified specific thermal conductivity z.",
            )
            if dim == 3:
                _require_keys(
                    state,
                    ["thermal_conductivity_y"],
                    "Each state needs a user-specified specific thermal conductivity y.",
                )
            if uses_convective:
                _require_keys(
                    state,
                    ["convection_heat_transfer_coef"],
                    "Convective BCs require a user-specified convection heat transfer coefficient.",
                )
            if uses_radiative:
                _require_keys(
                    state,
                    ["emissivity"],
                    "Radiative BCs require a user-specified emissivity.",
                )

        if not has_a_valid_state:
            raise ValidationError(
                "Material without any valid state (solid, powder, or liquid)."
            )
        if uses_convective:
            _require_keys(
                material,
                ["convection_temperature_infty"],
                "Convective BCs require setting 'convection_temperature_infty'.",
            )
        if uses_radiative:
            _require_keys(
                material,
                ["radiation_temperature_infty"],
                "Radiative BCs require setting 'radiation_temperature_infty'.",
            )

    # Tree: memory_space
    memory_space = _get(database, "memory_space", _OPTIONAL_MISSING)
    if memory_space is not _OPTIONAL_MISSING:
        memory_space = _as_string(memory_space, "memory_space")
        if memory_space not in {"device", "host"}:
            raise ValidationError(
                f"Method type for memory space, '{memory_space}', is not recognized. "
                "Valid options are: 'host' and 'device'"
            )

    # Tree: post_processor
    post_processor = _require_object(database, "post_processor")
    if "filename_prefix" not in post_processor:
        raise ValidationError(
            "The filename prefix for the postprocessor must be specified."
        )

    # Tree: refinement
    if "refinement" not in database:
        raise ValidationError("A refinement section of the input file must exist.")

    # Tree: sources
    sources = _require_object(database, "sources")
    n_beams = _as_unsigned(_get(database, "sources.n_beams"), "sources.n_beams")
    for beam_index in range(n_beams):
        beam_name = f"beam_{beam_index}"
        beam = _require_object(sources, beam_name)
        beam_type = _as_string(
            _get(beam, "type"), f"sources.{beam_name}.type"
        )
        if beam_type.lower() not in {"goldak", "electron_beam", "cube", "gaussian"}:
            raise ValidationError(
                f"Beam type, '{beam_type}', is not recognized. Valid options are: "
                "'goldak', 'electron_beam', 'cube', and 'gaussian'."
            )
        if "scan_path_file" not in beam:
            raise ValidationError("A scan path file for each beam must be given.")
        file_format = _as_string(
            _get(beam, "scan_path_file_format"),
            f"sources.{beam_name}.scan_path_file_format",
        )
        if file_format.lower() not in {"segment", "event_series"}:
            raise ValidationError(
                f"Scan path file format, '{file_format}', is not recognized. "
                "Valid options are: 'segment' and 'event_series'."
            )
        depth = _as_number(_get(beam, "depth"), f"sources.{beam_name}.depth")
        if depth < 0.0:
            raise ValidationError("Heat source depth must be non-negative.")
        absorption_efficiency = _as_number(
            _get(beam, "absorption_efficiency"),
            f"sources.{beam_name}.absorption_efficiency",
        )
        if not 0.0 <= absorption_efficiency <= 1.0:
            raise ValidationError(
                "Heat source absorption efficiency must be between 0 and 1."
            )

    # Tree: time_stepping
    time_stepping_method = _as_string(
        _get(database, "time_stepping.method"), "time_stepping.method"
    )
    if time_stepping_method.lower() not in {
        "forward_euler",
        "rk_third_order",
        "rk_fourth_order",
    }:
        raise ValidationError(
            f"Time stepping method, '{time_stepping_method}', is not recognized. "
            "Valid options are: 'forward_euler', 'rk_third_order', and "
            "'rk_fourth_order'."
        )

    scan_path_for_duration = _as_bool(
        _get(database, "time.scan_path_for_duration", False),
        "time.scan_path_for_duration",
    )
    if scan_path_for_duration:
        duration = _as_number(
            _get(database, "time_stepping.duration"), "time_stepping.duration"
        )
        if duration < 0.0:
            raise ValidationError("Time stepping duration must be non-negative.")

    time_step = _as_number(
        _get(database, "time_stepping.time_step"), "time_stepping.time_step"
    )
    if time_step < 0.0:
        raise ValidationError("Time step must be non-negative.")

    # Tree: experiment
    experiment = _get(database, "experiment", _OPTIONAL_MISSING)
    if experiment is not _OPTIONAL_MISSING and isinstance(experiment, dict):
        experiment_active = _as_bool(
            _get(experiment, "read_in_experimental_data", False),
            "experiment.read_in_experimental_data",
        )
        if experiment_active:
            if "file" not in experiment:
                raise ValidationError(
                    "If reading experimental data, a file must be given."
                )
            if "last_frame" not in experiment:
                raise ValidationError(
                    "If reading experimental data, a last frame index must be given."
                )
            experiment_format = _as_string(
                _get(experiment, "format"), "experiment.format"
            )
            if experiment_format.lower() not in {"point_cloud", "ray"}:
                raise ValidationError(
                    "Experiment format must be 'point_cloud' or 'ray'."
                )
            first_frame_index = _as_unsigned(
                _get(experiment, "first_frame", 0), "experiment.first_frame"
            )
            last_frame_index = _as_unsigned(
                _get(experiment, "last_frame"), "experiment.last_frame"
            )
            if last_frame_index < first_frame_index:
                raise ValidationError(
                    "When reading experimental data, the last frame index cannot "
                    "be lower than the first frame index."
                )
            first_camera_id = _as_unsigned(
                _get(experiment, "first_camera_id"), "experiment.first_camera_id"
            )
            last_camera_id = _as_unsigned(
                _get(experiment, "last_camera_id"), "experiment.last_camera_id"
            )
            if last_camera_id < first_camera_id:
                raise ValidationError(
                    "When reading experimental data, the last camera id cannot "
                    "be lower than the first camera id."
                )
            if "log_filename" not in experiment:
                raise ValidationError(
                    "If reading experimental data, a log filename must be given."
                )

    # Tree: data_assimilation
    convergence_tolerance = _get(
        database, "data_assimilation.convergence_tolerance", _OPTIONAL_MISSING
    )
    if convergence_tolerance is not _OPTIONAL_MISSING:
        if _as_number(
            convergence_tolerance, "data_assimilation.convergence_tolerance"
        ) < 0.0:
            raise ValidationError(
                "The data assimilation convergene tolerance must be non-negative."
            )

    localization_cutoff_function = _as_string(
        _get(database, "data_assimilation.localization_cutoff_function", "none"),
        "data_assimilation.localization_cutoff_function",
    )
    if localization_cutoff_function.lower() not in {
        "gaspari_cohn",
        "step_function",
        "none",
    }:
        raise ValidationError(
            "Unknown localization cutoff function. Valid options are "
            "'gaspari_cohn', 'step_function', and 'none'."
        )

    # Tree: units
    mesh_unit = _get(database, "units.mesh", _OPTIONAL_MISSING)
    if mesh_unit is not _OPTIONAL_MISSING:
        mesh_unit = _as_string(mesh_unit, "units.mesh")
        if mesh_unit.lower() not in {"millimeter", "centimeter", "inch", "meter"}:
            raise ValidationError(
                "Unknown unit associated with the mesh. Valid options are "
                "`millimeter`, `centimeter`, `inch`, and `meter`"
            )

    heat_source_power_unit = _get(
        database, "units.heat_source.power", _OPTIONAL_MISSING
    )
    if heat_source_power_unit is not _OPTIONAL_MISSING:
        heat_source_power_unit = _as_string(
            heat_source_power_unit, "units.heat_source.power"
        )
        if heat_source_power_unit.lower() not in {"milliwatt", "watt"}:
            raise ValidationError(
                "Unknown unit associated with the power of the heat source. "
                "Valid options are `milliwatt`, and `watt`"
            )

    heat_source_velocity_unit = _get(
        database, "units.heat_source.velocity", _OPTIONAL_MISSING
    )
    if heat_source_velocity_unit is not _OPTIONAL_MISSING:
        heat_source_velocity_unit = _as_string(
            heat_source_velocity_unit, "units.heat_source.velocity"
        )
        if heat_source_velocity_unit.lower() not in {
            "millimeter/second",
            "centimeter/second",
            "meter/second",
        }:
            raise ValidationError(
                "Unknown unit associated with the velocity of the heat source. "
                "Valid options are `millimeter/second`, `centimeter/second`, "
                "and `meter/second`"
            )

    heat_source_dimension_unit = _get(
        database, "units.heat_source.dimension", _OPTIONAL_MISSING
    )
    if heat_source_dimension_unit is not _OPTIONAL_MISSING:
        heat_source_dimension_unit = _as_string(
            heat_source_dimension_unit, "units.heat_source.dimension"
        )
        if heat_source_dimension_unit.lower() not in {
            "millimeter",
            "centimeter",
            "inch",
            "meter",
        }:
            raise ValidationError(
                "Unknown unit associated with the dimension of the heat source. "
                "Valid options are `millimeter`, `centimeter`, `inch`, and `meter`"
            )

    heat_source_scan_path_unit = _get(
        database, "units.heat_source.scan_path", _OPTIONAL_MISSING
    )
    if heat_source_scan_path_unit is not _OPTIONAL_MISSING:
        heat_source_scan_path_unit = _as_string(
            heat_source_scan_path_unit, "units.heat_source.scan_path"
        )
        if heat_source_scan_path_unit.lower() not in {
            "millimeter",
            "centimeter",
            "inch",
            "meter",
        }:
            raise ValidationError(
                "Unknown unit associated with the scan path. Valid options are "
                "`millimeter`, `centimeter`, `inch`, and `meter`"
            )


def _reject_nonstandard_json(value: str) -> None:
    raise ValueError(f"Non-standard JSON constant {value!r} is not allowed")


def _load_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as input_file:
            database = json.load(input_file, parse_constant=_reject_nonstandard_json)
    except OSError as error:
        raise ValidationError(f"Cannot read '{path}': {error}") from error
    except json.JSONDecodeError as error:
        raise ValidationError(
            f"Invalid JSON in '{path}' at line {error.lineno}, column "
            f"{error.colno}: {error.msg}"
        ) from error
    except ValueError as error:
        raise ValidationError(f"Invalid JSON in '{path}': {error}") from error

    if not isinstance(database, dict):
        raise ValidationError("The top-level JSON value must be an object.")
    return database


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate an Adamantine JSON input file."
    )
    parser.add_argument("input_file", type=Path, help="path to the JSON input file")
    args = parser.parse_args(argv)

    try:
        database = _load_json(args.input_file)
        validate_input_database(database)
    except ValidationError as error:
        print(f"Invalid Adamantine input: {error}", file=sys.stderr)
        return 1

    print(f"Valid Adamantine input: {args.input_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
