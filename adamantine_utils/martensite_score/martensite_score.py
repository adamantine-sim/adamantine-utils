"""Score an Adamantine simulation by the time points spend in the precipitation band after becoming martensite.

A point becomes martensite when it cools through the martensite-start temperature (Ms) and stays martensite until it
is reheated above the re-austenitization temperature. The per-point score is the number of time steps spent inside the
precipitation band while martensitic. The total score is the sum over all points. Temperatures are in Kelvin.
"""

import argparse
import os
from dataclasses import dataclass, replace

import numpy as np
import pyvista

from adamantine_utils.melt_pool_analysis.melt_pool_analysis import get_iteration_count

K_TO_C = 273.15


@dataclass(frozen=True)
class MaterialParameters:
    """Critical temperatures of a material (Kelvin)."""

    martensite_start: float  # Ms
    reaustenization: float  # Above this temperature martensite is destroyed
    precip_center: float  # Optimal precipitation temperature
    precip_width: float  # Full width of the band centered on precip_center

    def __post_init__(self):
        if not self.martensite_start < self.precip_min:
            raise ValueError("The precipitation band must lie above the martensite start temperature.")
        if not self.precip_max <= self.reaustenization:
            raise ValueError("The precipitation band must not exceed the re-austenitization temperature.")

    @property
    def precip_min(self):
        return self.precip_center - self.precip_width / 2.0

    @property
    def precip_max(self):
        return self.precip_center + self.precip_width / 2.0

    @classmethod
    def from_celsius(cls, martensite_start, reaustenization, precip_center, precip_width):
        return cls(martensite_start + K_TO_C, reaustenization + K_TO_C, precip_center + K_TO_C, precip_width)


MATERIALS = {
    "17-4PH": MaterialParameters.from_celsius(
        martensite_start=105.0, reaustenization=565.0, precip_center=522.5, precip_width=42.5
    ),
}


def _sort_order(points):
    """Order that sorts points lexicographically by (x, y, z)."""
    return np.lexsort(points.T[::-1])


def get_time_series(directory, base_filename, field_name="temperature"):
    """Return (iteration_numbers, array of shape (n_times, n_points)).

    Points are sorted by (x, y, z) so that they line up between time steps.
    """
    iteration_numbers = get_iteration_count(os.path.join(directory, base_filename))
    if not iteration_numbers:
        raise FileNotFoundError(f"No {base_filename}.*.pvtu files found in {directory}")

    values = []
    ref_points = None
    for iteration in iteration_numbers:
        dataset = pyvista.read(os.path.join(directory, f"{base_filename}.{iteration}.pvtu"))
        order = _sort_order(dataset.points)
        points = dataset.points[order]
        if ref_points is None:
            ref_points = points
        elif points.shape != ref_points.shape:
            raise ValueError(f"Point count mismatch at iteration {iteration}")
        elif not np.allclose(points, ref_points, atol=1e-10):
            print(f"WARNING: iteration {iteration} points do not match the first time step")
        values.append(dataset.point_data[field_name][order])

    return iteration_numbers, np.array(values)


def martensite_score_per_point(time_series, params):
    """Compute the score of each point from a (n_times, n_points) temperature array."""
    martensite = np.zeros(time_series.shape[1], dtype=bool)
    scores = np.zeros(time_series.shape[1])
    for previous, current in zip(time_series[:-1], time_series[1:]):
        # Martensite forms on cooling through Ms and is destroyed by reheating above the re-austenitization temperature
        martensite |= (previous > params.martensite_start) & (current < params.martensite_start)
        martensite &= current <= params.reaustenization
        scores += martensite & (current >= params.precip_min) & (current <= params.precip_max)
    return scores


def martensite_score(directory, base_filename, params, field_name="temperature", line_plots=False,
                     volume_plot=False):
    """Return the total martensite score of a simulation."""
    iteration_numbers, time_series = get_time_series(directory, base_filename, field_name)
    per_point_scores = martensite_score_per_point(time_series, params)
    if line_plots:
        plot_point_histories(iteration_numbers, time_series)
    if volume_plot:
        plot_score_on_mesh(directory, base_filename, iteration_numbers[0], per_point_scores)
    return float(per_point_scores.sum())


def plot_point_histories(iteration_numbers, time_series, n_points=5):
    import matplotlib.pyplot as plt

    for i in range(min(n_points, time_series.shape[1])):
        plt.figure()
        plt.plot(iteration_numbers, time_series[:, i], "b.-")
        plt.title(f"Point {i} history")
        plt.xlabel("Iteration")
        plt.ylabel("Temperature (K)")
        plt.grid(True)
    plt.show()


def plot_score_on_mesh(directory, base_filename, first_iteration, per_point_scores):
    dataset = pyvista.read(os.path.join(directory, f"{base_filename}.{first_iteration}.pvtu"))
    # per_point_scores follows the sorted point order used by get_time_series
    scores = np.empty_like(per_point_scores)
    scores[_sort_order(dataset.points)] = per_point_scores
    dataset.point_data["score"] = scores
    dataset.set_active_scalars("score")
    plotter = pyvista.Plotter()
    plotter.add_mesh(dataset.clip("y"), show_edges=True)
    plotter.show()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data-directory", required=True, help="location of the VTK files")
    parser.add_argument("-f", "--base-filename", required=True, help="base filename of the VTK files")
    parser.add_argument("-v", "--variable", default="temperature", help="temperature variable name in the files")
    parser.add_argument("--material", choices=sorted(MATERIALS), help="material preset")
    parser.add_argument("--ms", type=float, help="martensite start temperature (C)")
    parser.add_argument("--reaustenization", type=float, help="re-austenitization temperature (C)")
    parser.add_argument("--precip-center", type=float, help="optimal precipitation temperature (C)")
    parser.add_argument("--precip-width", type=float, help="full width of the precipitation band (C)")
    parser.add_argument("--line-plots", action="store_true", help="plot some single-point temperature histories")
    parser.add_argument("--volume-plot", action="store_true", help="plot a 3D cross-section of the per-point score")
    args = parser.parse_args(argv)

    overrides = {"martensite_start": args.ms, "reaustenization": args.reaustenization,
                 "precip_center": args.precip_center, "precip_width": args.precip_width}
    overrides = {k: v for k, v in overrides.items() if v is not None}
    if args.material:
        # Overrides are given in Celsius, except the width, which is a difference
        params = replace(MATERIALS[args.material],
                         **{k: v if k == "precip_width" else v + K_TO_C for k, v in overrides.items()})
    elif len(overrides) == 4:
        params = MaterialParameters.from_celsius(**overrides)
    else:
        parser.error("Give --material or all of --ms, --reaustenization, --precip-center and --precip-width")

    score = martensite_score(args.data_directory, args.base_filename, params, args.variable,
                             args.line_plots, args.volume_plot)
    print("Score:", score)
    return score


if __name__ == "__main__":
    main()
