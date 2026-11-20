import os

import numpy as np
import pytest

from adamantine_utils.martensite_score.martensite_score import (
    MATERIALS, MaterialParameters, main, martensite_score, martensite_score_per_point)

# Ms=100, T_A=500, band [300, 400] (Kelvin)
P = MaterialParameters(100.0, 500.0, 350.0, 100.0)


def score(*columns):
    return martensite_score_per_point(np.array(columns, dtype=float).T, P)


def test_cooling_then_band():
    # heat above Ms, cool below Ms, reheat into band twice
    assert score([20, 200, 50, 350, 350, 20])[0] == 2


def test_never_above_ms():
    assert score([20, 50, 350, 350])[0] == 0


def test_band_before_martensite_not_counted():
    assert score([20, 350, 350, 50, 350])[0] == 1


def test_reaustenization_resets():
    assert score([20, 200, 50, 600, 350, 350])[0] == 0


def test_validation():
    with pytest.raises(ValueError):
        MaterialParameters(100.0, 500.0, 120.0, 50.0)
    with pytest.raises(ValueError):
        MaterialParameters(100.0, 500.0, 480.0, 100.0)


def original_17_4(time_series):
    T_Ms, T_ppt, T_band, T_A = 378.15, 795.65, 42.5, 838.15
    T, N = time_series.shape
    seen = np.zeros(N, dtype=bool)
    on = np.zeros(N, dtype=bool)
    scores = np.zeros(N)
    prev = time_series[0]
    for n in range(1, T):
        cur = time_series[n]
        seen |= cur > T_Ms
        on |= (prev > T_Ms) & (cur < T_Ms) & seen
        on &= ~(cur > T_A)
        in_band = (cur >= T_ppt - T_band / 2) & (cur <= T_ppt + T_band / 2)
        scores += on * in_band
        prev = cur
    return scores


def test_matches_original_17_4():
    rng = np.random.default_rng(0)
    series = rng.uniform(250, 900, size=(200, 50))
    series[0] = 293.15  # simulations start below Ms
    np.testing.assert_array_equal(martensite_score_per_point(series, MATERIALS["17-4PH"]), original_17_4(series))


def test_on_vtu_files():
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_vtu_files")
    result = martensite_score(path, "solution", MATERIALS["17-4PH"])
    assert np.isfinite(result) and result >= 0


def write_run(path, temperatures, shuffle_last=False):
    """Write a one-piece-per-step run with points at x=0,1,... and the given (n_times, n_points) temperatures."""
    import pyvista

    n_points = temperatures.shape[1]
    for step, values in enumerate(temperatures):
        points = np.zeros((n_points, 3))
        points[:, 0] = np.arange(n_points)
        values = values.astype(np.float32)
        if shuffle_last and step % 2:  # store points in reversed order on odd steps
            points, values = points[::-1], values[::-1]
        grid = pyvista.PolyData(points.copy()).cast_to_unstructured_grid()
        grid.point_data["temperature"] = values
        grid.save(str(path / f"run.{step}.0.vtu"))
        (path / f"run.{step}.pvtu").write_text(
            '<?xml version="1.0"?>\n<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
            '<PUnstructuredGrid GhostLevel="0"><PPointData><PDataArray type="Float32" Name="temperature"/>'
            '</PPointData><PPoints><PDataArray type="Float32" NumberOfComponents="3"/></PPoints>'
            f'<Piece Source="run.{step}.0.vtu"/></PUnstructuredGrid></VTKFile>\n')


def test_known_score_from_files(tmp_path):
    # Point 0 becomes martensite and sits in the band for 2 steps; point 1 never leaves ambient
    temperatures = np.array([[20, 20], [200, 20], [50, 20], [350, 20], [350, 20]], dtype=float)
    write_run(tmp_path, temperatures, shuffle_last=True)
    assert martensite_score(str(tmp_path), "run", P) == 2


def test_missing_files(tmp_path):
    with pytest.raises(FileNotFoundError):
        martensite_score(str(tmp_path), "run", P)


def test_cli_preset_and_override(tmp_path, capsys):
    # Peak at 400 C, cool to 20 C, then 300 C: in the 17-4PH band only if its center is moved to 300 C
    temperatures = np.array([[20], [400], [20], [300], [300]], dtype=float) + 273.15
    write_run(tmp_path, temperatures)
    args = ["-d", str(tmp_path), "-f", "run"]
    assert main(args + ["--material", "17-4PH"]) == 0
    assert main(args + ["--material", "17-4PH", "--precip-center", "300"]) == 2
    assert main(args + ["--ms", "105", "--reaustenization", "565", "--precip-center", "300",
                        "--precip-width", "42.5"]) == 2
    assert "Score: 2.0" in capsys.readouterr().out


def test_cli_requires_all_parameters(tmp_path):
    with pytest.raises(SystemExit):
        main(["-d", str(tmp_path), "-f", "run", "--ms", "105"])
