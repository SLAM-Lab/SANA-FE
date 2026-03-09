#!/usr/bin/env python3
"""
    python3 scripts/test_performance_viz.py
"""
from pathlib import Path
import sys
import importlib.util

REPO = Path(__file__).resolve().parent.parent


def _bootstrap_sanafe_viz():
    import types
    for pkg in ("sanafe", "sanafe.data", "sanafe.viz"):
        if pkg not in sys.modules:
            mod = types.ModuleType(pkg)
            mod.__path__ = []
            sys.modules[pkg] = mod

    def load(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    load("sanafe.data.traces", REPO / "sanafe" / "data" / "traces.py")
    load("sanafe.viz.styles", REPO / "sanafe" / "viz" / "styles.py")
    load("sanafe.viz.performance", REPO / "sanafe" / "viz" / "performance.py")

    from sanafe.viz.performance import (
        energy_breakdown_plot, throughput_plot,
        latency_histogram, latency_comparison,
    )
    from sanafe.data.traces import TraceData
    return (TraceData, energy_breakdown_plot, throughput_plot,
            latency_histogram, latency_comparison)


def main():
    import matplotlib.pyplot as plt

    (TraceData, energy_breakdown_plot, throughput_plot,
     latency_histogram, latency_comparison) = _bootstrap_sanafe_viz()

    out = REPO / "tmp_trace_test"
    out.mkdir(exist_ok=True)

    # --- Load the tutorial traces that already exist on disk ---
    perf_csv = REPO / "tutorial" / "perf.csv"
    msg_csv = REPO / "tutorial" / "messages.csv"
    tr = TraceData.from_files(perf_csv=perf_csv, message_csv=msg_csv)
    print("TraceData:", tr)

    # 1. Energy breakdown — stacked area (default)
    fig, ax = energy_breakdown_plot(tr, title="Energy Breakdown (stacked area)")
    fig.savefig(out / "energy_stacked_area.png", dpi=100)
    print("Saved", out / "energy_stacked_area.png")

    # 2. Energy breakdown — stacked bar
    fig, ax = energy_breakdown_plot(tr, mode="stacked_bar", title="Energy Breakdown (stacked bar)")
    fig.savefig(out / "energy_stacked_bar.png", dpi=100)
    print("Saved", out / "energy_stacked_bar.png")

    # 3. Energy breakdown — normalized percentage
    fig, ax = energy_breakdown_plot(tr, mode="stacked_bar", normalize=True,
                                    title="Energy Breakdown (normalized)")
    fig.savefig(out / "energy_normalized.png", dpi=100)
    print("Saved", out / "energy_normalized.png")

    # 4. Throughput
    fig, ax = throughput_plot(tr, title="Throughput")
    fig.savefig(out / "throughput.png", dpi=100)
    print("Saved", out / "throughput.png")

    # 5. Throughput — dual axis with energy
    fig, ax = throughput_plot(tr, metrics=["fired", "total_energy"],
                              secondary_y=["total_energy"],
                              title="Fired vs Total Energy")
    fig.savefig(out / "throughput_dual.png", dpi=100)
    print("Saved", out / "throughput_dual.png")

    # 6. Latency histogram — generation delay
    fig, ax = latency_histogram(tr, metric="generation_delay",
                                title="Generation Delay")
    fig.savefig(out / "latency_generation.png", dpi=100)
    print("Saved", out / "latency_generation.png")

    # 7. Latency histogram — processing delay
    fig, ax = latency_histogram(tr, metric="processing_delay",
                                title="Processing Delay")
    fig.savefig(out / "latency_processing.png", dpi=100)
    print("Saved", out / "latency_processing.png")

    # 8. Latency comparison overlay
    fig, ax = latency_comparison(tr, metrics=["generation_delay", "processing_delay"])
    fig.savefig(out / "latency_comparison.png", dpi=100)
    print("Saved", out / "latency_comparison.png")

    plt.close("all")
    print("Done. Check files in tmp_trace_test/")


if __name__ == "__main__":
    main()
