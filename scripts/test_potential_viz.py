#!/usr/bin/env python3
"""
    python3 scripts/test_potential_viz.py
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
    load("sanafe.viz.potential", REPO / "sanafe" / "viz" / "potential.py")

    from sanafe.viz.potential import potential_plot, potential_heatmap, potential_subplots
    from sanafe.data.traces import TraceData
    return TraceData, potential_plot, potential_heatmap, potential_subplots


def main():
    import pandas as pd
    import numpy as np

    TraceData, potential_plot, potential_heatmap, potential_subplots = _bootstrap_sanafe_viz()

    out = REPO / "tmp_trace_test"
    out.mkdir(exist_ok=True)

    pot_csv = out / "potentials.csv"
    neuron_ids = ["in.0", "in.1", "out.0"]
    pd.DataFrame({
        "in.0": [0.0, 0.1, 0.2, 0.3, 0.25, 0.1],
        "in.1": [0.0, 0.0, 0.5, 0.0, 0.2, 0.0],
        "out.0": [0.0, 0.0, 0.0, 0.9, 0.1, 0.0],
    }).to_csv(pot_csv, index=False)

    tr = TraceData.from_files(potential_csv=pot_csv)
    print("TraceData:", tr)

    # Timeseries plot
    fig1, ax1 = potential_plot(
        tr,
        neuron_ids=neuron_ids,
        title="Potential test (from CSV)",
        show_threshold=0.5,
    )
    fig1.savefig(out / "potential_plot.png", dpi=100)
    print("Saved", out / "potential_plot.png")

    # Heatmap
    fig2, ax2 = potential_heatmap(
        tr,
        neuron_ids=neuron_ids,
        title="Potential heatmap (from CSV)",
    )
    fig2.savefig(out / "potential_heatmap.png", dpi=100)
    print("Saved", out / "potential_heatmap.png")

    # Subpotentials
    fig3, axes = potential_subplots(tr, neuron_ids=neuron_ids, ncols=1)
    fig3.savefig(out / "potential_subplots.png", dpi=100)
    print("Saved", out / "potential_subplots.png")

    # Heatmap from array 
    arr = np.random.rand(20, 3).astype(np.float64) * 0.8
    fig4, ax4 = potential_heatmap(arr, neuron_ids=["A", "B", "C"], title="From numpy array")
    fig4.savefig(out / "potential_heatmap_array.png", dpi=100)
    print("Saved", out / "potential_heatmap_array.png")

    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except Exception:
        pass
    print("Done. Check files in tmp_trace_test/")


if __name__ == "__main__":
    main()
