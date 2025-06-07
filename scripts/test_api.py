import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
import sanafe

arch = sanafe.load_arch("arch/example.yaml")
chip = sanafe.SpikingChip(arch)

net = sanafe.load_net("snn/example.net", arch, use_netlist_format=True)
net.save("out", use_netlist_format=True)
chip.load(net)

result = chip.sim(10, message_trace=True)

import yaml
print(yaml.dump(result))
