import importlib.util,os,sys,numpy as np,pandas as pd
P=os.getcwd(); out=P+'/runs/neuroscale/data/'
spec=importlib.util.spec_from_file_location("ns",P+"/scripts/ongoing/neuroscale.py"); ns=importlib.util.module_from_spec(spec); spec.loader.exec_module(ns)
arch=ns.make_arch_with_protocol(ns.load_arch_dict(P+'/arch/neuroscale3.yaml'),'neuroscale')
tiles,cm=ns.compute_tile_order(arch)
m=pd.read_csv(out+'messages.csv'); m=m[m.mid>=0]
pc=lambda hw: cm[int(str(hw).split('.')[0])*4+int(str(hw).split('.')[1])]
M=np.zeros((100,64,64))
for s,a,b,sp in zip(m.timestep,m.src_hw.map(pc),m.dest_hw.map(pc),m.spikes): M[s-1,a,b]+=1
np.save(out+'M.npy',M); print(M.sum(), M.sum(axis=(0,1)).reshape(16,4).sum(1))
