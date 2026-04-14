from llama_cpp import Llama
import json
import yaml
import os
import time
import sys

MODEL_PATH = "/home/usr1/jboyle/neuro/llama.cpp/models/gemma-4-26B-A4B-it-UD-Q4_K_M.gguf"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
# Setup paths
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

print("Initializing LLM...")
llm = Llama(
    model_path=MODEL_PATH,
    n_threads=12,
    n_ctx=4096,
    verbose=False,
)

SYSTEM_PROMPT = (
    "You are an optimization agent tuning integer parameter N (256≥N≥1) to minimize latency with the smallest value of N. "
    "Try to also minimize the number of steps required to get the optimal solution. "
    "If you find a plateau, continue to explore in the vicinity until you are confident. "
    "The simulation is deterministic - identical N always produces identical latency. Never retest a value. "
    "You will receive a table of all previously tested values and their latencies as relative differences. "
    "Respond with a JSON object containing:\n"
    '  "thought": your reasoning,\n'
    '  "next_n": the next integer value of N to test,\n'
    '  "status": "CONTINUE" or "STOP"\n'
    "Set STOP only when you are confident further testing will not find a better N."
)

def format_history(history):
    if not history:
        return "No results yet. Suggest an initial value of N to test."

    best_lat = min(h["latency"] for h in history)
    rows = sorted(history, key=lambda h: h["latency"])

    lines = ["N | delta"]
    for h in rows:
        delta = (h["latency"] - best_lat) / best_lat
        lines.append(f"{h['n']} | +{delta:.2%}")

    lines.append(f"\nBest: N={rows[0]['n']}")
    lines.append(f"Tested: {len(history)}")
    return "\n".join(lines)


def build_context(history):
    """Build a compact context string for the LLM."""
    if not history:
        return "No results yet. Suggest an initial value of N to test."

    formatted = format_history(history)
    return formatted


def create_modified_arch_file(original_arch_path, max_parallel_accesses, run_dir):
    with open(original_arch_path, "r") as f:
        arch_data = yaml.safe_load(f)

    for tile in arch_data["architecture"]["tile"]:
        for core in tile["core"]:
            for synapse in core["synapse"]:
                if synapse["name"] == "loihi_conv_synapse":
                    synapse["attributes"]["plugin"] = (
                        "/home/usr1/jboyle/neuro/sana-fe/plugins/libloihi_synapse.so"
                    )
                    synapse["attributes"]["model"] = "loihi"
                    synapse["attributes"]["max_parallel_accesses"] = max_parallel_accesses
                    synapse["attributes"].pop("latency_process_spike", None)
                    synapse["attributes"].pop("energy_process_spike", None)

    new_arch_filename = f"loihi_parallel_{max_parallel_accesses}.yaml"
    new_arch_path = os.path.join(run_dir, new_arch_filename)
    with open(new_arch_path, "w") as f:
        yaml.dump(arch_data, f, default_flow_style=False)
    return new_arch_path


def run_simulation(arch_path, network_path, timesteps):
    arch = sanafe.load_arch(arch_path)
    net = sanafe.load_net(network_path, arch, use_netlist_format=True)
    chip = sanafe.SpikingChip(arch)
    chip.load(net)

    start = time.perf_counter()
    result = chip.sim(
        timesteps,
        timing_model="simple",
        processing_threads=8,
        scheduler_threads=1,
    )
    chip.reset()
    wall_time = time.perf_counter() - start
    return result["sim_time"], wall_time


def query_llm(llm, history):
    """Send compressed context to LLM, get back a single next N."""
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": build_context(history)},
    ]

    output = llm.create_chat_completion(
        messages=messages,
        response_format={
            "type": "json_object",
            "schema": {
                "type": "object",
                "properties": {
                    "thought": {"type": "string"},
                    "next_n": {"type": "integer"},
                    "status": {
                        "type": "string",
                        "enum": ["CONTINUE", "STOP"],
                    },
                },
                "required": ["thought", "next_n", "status"],
            },
        },
    )

    res = json.loads(output["choices"][0]["message"]["content"])
    return res


# --- Paths ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))
RUN_DIR = os.path.join(PROJECT_DIR, "runs", "llm")

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", "loihi.yaml")
NETWORK_PATH = os.path.join(RUN_DIR, "dvs_gesture_32x32.net.tagged")

TIMESTEPS = 128
MAX_ITERATIONS = 50

# --- Main loop ---
history = []
tested_ns = set()
retries = 0
MAX_RETRIES = 3

for iteration in range(MAX_ITERATIONS):
    res = query_llm(llm, history)
    print(f"\n--- Iteration {iteration + 1} ---")
    print(f"Thought: {res['thought']}")
    print(f"Next N: {res['next_n']}")
    print(f"Status: {res['status']}")

    if res["status"] == "STOP":
        if history:
            best = min(history, key=lambda h: h["latency"])
            print(f"\nAgent stopped. Best N={best['n']} latency={best['latency']:.6e}")
        break

    n = res["next_n"]
    if n in tested_ns or n < 1:
        retries += 1
        print(f"Warning: N={n} already tested or invalid (retry {retries}/{MAX_RETRIES})")
        if retries >= MAX_RETRIES:
            print("Max retries on duplicate values. Stopping.")
            break
        continue

    retries = 0
    print(f"  Simulating N={n}...", end=" ", flush=True)
    arch_path = create_modified_arch_file(ARCH_PATH, n, RUN_DIR)
    latency, wall = run_simulation(arch_path, NETWORK_PATH, TIMESTEPS)
    history.append({"n": n, "latency": latency})
    tested_ns.add(n)
    print(f"latency={latency:.6e} (wall={wall:.2f}s)")

best = min(history, key=lambda h: h["latency"])
print(f"\n  Final result: Best N={best['n']} with latency={best['latency']:.6e}")
print(f"   Tested {len(history)} unique values of N")
