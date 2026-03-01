from pathlib import Path
import sys
import importlib.util

REPO = Path(__file__).resolve().parent.parent


def _bootstrap_sanafe_viz():
    """Load sanafe.data.traces and sanafe.viz.raster without running sanafe/__init__.py."""
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
    load("sanafe.viz.raster", REPO / "sanafe" / "viz" / "raster.py")

    from sanafe.viz.raster import raster_plot, raster_plot_matrix
    from sanafe.data.traces import TraceData
    return TraceData, raster_plot, raster_plot_matrix


def main():
    import pandas as pd

    TraceData, raster_plot, raster_plot_matrix = _bootstrap_sanafe_viz()

    out = REPO / "tmp_trace_test"
    out.mkdir(exist_ok=True)

    spikes_csv = out / "spikes.csv"
    pd.DataFrame([
        {"neuron": "in.0", "timestep": 0},
        {"neuron": "in.1", "timestep": 0},
        {"neuron": "in.0", "timestep": 1},
        {"neuron": "in.1", "timestep": 2},
        {"neuron": "out.0", "timestep": 2},
        {"neuron": "out.0", "timestep": 3},
        {"neuron": "out.1", "timestep": 4},
    ]).to_csv(spikes_csv, index=False)

    tr = TraceData.from_files(spike_csv=spikes_csv)
    print("TraceData:", tr)
    print("Groups:", tr.get_neuron_groups())

    fig1, ax1 = raster_plot(tr, title="Raster (from CSV, all groups)")
    fig1.savefig(out / "raster_plot.png", dpi=100)
    print("Saved", out / "raster_plot.png")

    fig2, ax2 = raster_plot(
        tr,
        groups=["in", "out"],
        time_range=(0, 4),
        title="Raster (groups=in,out; time 0–4)",
    )
    fig2.savefig(out / "raster_plot_filtered.png", dpi=100)
    print("Saved", out / "raster_plot_filtered.png")

    matrix, neuron_ids = tr.spikes_to_matrix()
    fig3, ax3 = raster_plot_matrix(
        matrix,
        neuron_ids=neuron_ids,
        title="Raster from matrix (spikes_to_matrix)",
    )
    fig3.savefig(out / "raster_plot_matrix.png", dpi=100)
    print("Saved", out / "raster_plot_matrix.png")

    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except Exception:
        pass
    print("Done. Check files in tmp_trace_test/")


if __name__ == "__main__":
    main()
