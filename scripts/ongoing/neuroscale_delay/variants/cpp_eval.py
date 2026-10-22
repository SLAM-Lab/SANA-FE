import importlib.util,os,sys,yaml,copy,numpy as np,pandas as pd
P=os.getcwd(); sys.path[:0]=[P+'/build_py',P]
sys.path.insert(0,'/tmp/claude-0/-opt-project-Neuromorphic-Architectures/2267edbb-a264-48d2-a1af-cbca8e496779/scratchpad/an')
os.environ['NEUROSCALE_LOG_DIR']=P+'/runs/neuroscale/cpp'
spec=importlib.util.spec_from_file_location("ns",P+"/scripts/ongoing/neuroscale.py"); ns=importlib.util.module_from_spec(spec); spec.loader.exec_module(ns)
import harness as H
base=ns.load_arch_dict(P+'/arch/neuroscale3.yaml')
def set_out(node,v):
    if isinstance(node,dict):
        for k in list(node):
            if k=='latency_message_out': node[k]=v
            else: set_out(node[k],v)
    elif isinstance(node,list):
        for x in node: set_out(x,v)
def run(name,extra,msg_out=None):
    d=copy.deepcopy(base); d['architecture']['attributes'].update(extra)
    if msg_out is not None: set_out(d,msg_out)
    d['architecture']['attributes']['sync_protocol']='neuroscale'
    os.makedirs(ns.LOG_DIR,exist_ok=True); path=ns.LOG_DIR+'/_a.yaml'; yaml.safe_dump(d,open(path,'w'),sort_keys=False)
    arch=ns.sanafe.load_arch(path); tiles,cm=ns.compute_tile_order(arch)
    net,_=ns.build_network_from_json(arch,ns.USE_JSON_NETWORK); chip=ns.sanafe.SpikingChip(arch); chip.load(net)
    mp=ns.LOG_DIR+'/m.csv'; chip.sim(100,timing_model='simple',message_trace=mp,spike_trace=False)
    m=pd.read_csv(mp); rows=ns.extract_neuroscale_timing(m,cm)
    T=np.zeros((64,100))
    for c,s,t in rows: T[c,s-1]=t
    T=T-T[:,[0]]
    a=H.score(T,1,100); b=H.score(T,50,100)
    print('%-42s MAPE %5.2f%% (test %5.2f%%) pop-finish err %5.2f%% spread %4.0fus (actsim %.0f)'%(name,a[0],b[0],a[1],a[2],a[3]))
    return T
if __name__=='__main__':
    run('C++ baseline (YAML as-is)',{})
    run('+ axon scan',{'neuroscale_axon_scan_latency':3.8e-9})
    run('+ tile port + sync pkts',{'neuroscale_axon_scan_latency':3.8e-9,'neuroscale_packet_latency':3.17e-9,'neuroscale_sync_packets':True},msg_out=3.17e-9)
    for Q in (40,50,66):
        run('+ stall Q=%d'%Q,{'neuroscale_axon_scan_latency':3.8e-9,'neuroscale_packet_latency':3.17e-9,'neuroscale_sync_packets':True,'neuroscale_stall_queue_depth':Q,'neuroscale_stall_drain_latency':11.2e-9},msg_out=3.17e-9)
