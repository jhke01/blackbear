#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def load_csv(csv_path):
    time = []
    avg_cluster_radius = []
    monomer_concentration = []
    total_cluster_density = []

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            time.append(float(row["time"]))
            avg_cluster_radius.append(float(row["avg_cluster_radius"]))
            monomer_concentration.append(float(row["monomer_concentration"]))
            total_cluster_density.append(float(row["total_cluster_density"]))

    return time, avg_cluster_radius, monomer_concentration, total_cluster_density


def main():
    parser = argparse.ArgumentParser(
        description="Plot cluster dynamics CSV output with total cluster density in #/m^3."
    )
    parser.add_argument(
        "csv_file",
        nargs="?",
        default="cluster_dynamics_50_interfacial_energy_out.csv",
        help="Path to the CSV file to plot.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output image path. Defaults to <csv_stem>_plots.png.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv_file).resolve()
    output_path = (
        Path(args.output).resolve()
        if args.output
        else csv_path.with_name(f"{csv_path.stem}_plots.png")
    )

    time, avg_cluster_radius, monomer_concentration, total_cluster_density = load_csv(csv_path)
    filtered = [
        (t, r, d)
        for t, r, d in zip(time, avg_cluster_radius, total_cluster_density)
        if t > 0.0
    ]
    time_plot = [item[0] for item in filtered]
    avg_cluster_radius_nm = [item[1] * 1e9 for item in filtered]
    total_cluster_density_number = [item[2] for item in filtered]

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    axes[0].plot(time_plot, total_cluster_density_number, color="black", linewidth=2)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlim(left=1e-2)
    axes[0].set_ylim(1e15, 1e27)
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Cluster number density (m$^{-3}$)")

    axes[1].plot(time_plot, avg_cluster_radius_nm, color="black", linewidth=2)
    axes[1].set_xscale("log")
    axes[1].set_xlim(left=1e-2)
    axes[1].set_ylim(0.0, 2.0)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Average cluster radius (nm)")
    axes[1].yaxis.set_minor_locator(AutoMinorLocator())

    for ax in axes:
        ax.tick_params(which="both", direction="in", top=True, right=True, labeltop=False)
        ax.tick_params(which="major", length=6)
        ax.tick_params(which="minor", length=3)
        ax.minorticks_on()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    print(f"Saved plot to {output_path}")


if __name__ == "__main__":
    main()
