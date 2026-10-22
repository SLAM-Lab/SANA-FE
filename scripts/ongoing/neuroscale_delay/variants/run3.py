import importlib; sm=importlib.import_module("sim2")
from sim2 import rep, P, n_fi, n_fo, NC, np
import harness
orig=sm.sync
def run(name,sy,port,scan=True):
    sm.sync=sy
    T,_=sm.sim2(P,scan,port); rep(name,T); sm.sync=orig
z=np.zeros(NC,int)
run('scan + per-core out(3ns) no sync',z,False)
run('scan + per-core out(3ns) + sync pkts only (3ns each)',n_fi+n_fo,False)
run('scan + tile port, no sync pkts',z,True)
run('scan + tile port + sync pkts',n_fi+n_fo,True)
