import json,numpy as np,pandas as pd
R='/tmp/claude-0/-opt-project-Neuromorphic-Architectures/2267edbb-a264-48d2-a1af-cbca8e496779/scratchpad/sana/SANA-FE/runs/neuroscale/data/'
NC,NS,TADV=64,100,6
act=pd.read_csv(R+'activity.csv')
def grid(col):
    return act.pivot(index='core',columns='step',values=col).reindex(index=range(NC),columns=range(1,NS+1)).fillna(0).to_numpy()
mo,mi,ev=grid('msgs_out'),grid('msgs_in'),grid('events_in')   # [core, i] i=SANA step-1
fan_in={int(k):v for k,v in json.load(open(R+'fan_in.json')).items()}
fan_out={int(k):v for k,v in json.load(open(R+'fan_out.json')).items()}
A=pd.read_csv(R+'actsim_T.csv',index_col=0).sort_index().to_numpy()
M=np.load(R+'M.npy')            # [i, src, dst]
Mt=M.reshape(NS,16,4,16,4).sum(axis=(2,))        # [i, srctile, dst tile, dcore]? -> (i,16,16,4)
n_fi=np.array([len(fan_in.get(c,())) for c in range(NC)]); n_fo=np.array([len(fan_out.get(c,())) for c in range(NC)])
FI=[np.array(sorted(fan_in.get(c,[])),int) for c in range(NC)]
FO=[np.array(sorted(fan_out.get(c,[])),int) for c in range(NC)]
dA=np.diff(A,axis=1)   # dA[:,j] <-> idx j+1
tile=np.arange(NC)//4
def tilesum(x): return np.repeat(x.reshape(16,4).sum(1),4)

def sim(p, variant='add'):
    """p: dict of params. Returns T[c, i] (i=0..99), plus stall flags."""
    T=np.zeros((NC,NS)); stall=np.zeros((NC,NS))
    for i in range(1,NS):
        prev=T[:,i-1]
        S=prev.copy()
        for c in range(NC):
            if len(FI[c]): S[c]=max(S[c],prev[FI[c]].max())
            if len(FO[c]) and i-TADV>=0: S[c]=max(S[c],T[FO[c],i-TADV].max())
        acc=p['t_ev']*ev[:,i-1]
        if variant=='add':
            L=p['F']+acc+p['t_in']*mi[:,i]+p['t_out']*mo[:,i]
        else:
            sync=n_fi+n_fo
            L=p['F']+acc+p['t_in']*mi[:,i]+p['t_pkt']*tilesum(mo[:,i]+sync)
        Tn=S+L
        if variant=='stall' or variant=='stall_tile':
            W=p['W']
            E=S+acc+W                 # sender emission start (neuron phase)
            U=S+acc+W                 # dest input opens (same instant for lockstep cores)
            # per sender tile x, dest core d: messages from tile x's cores to d
            for x in range(16):
                for d in np.nonzero(Mt[i,x].sum(0)>0)[0] if False else range(1):
                    pass
            # vectorised: Mt[i] shape (16 srctile,16 dsttile,4 dcore) -- recompute properly
            Ms=M[i].reshape(16,4,64).sum(1)   # [srctile, dst core]
            tin=p['t_in']; Q=p['Q']
            Et=E.reshape(16,4).max(1)     # tile emission start (HOL: tile moves together)
            for x in range(16):
                for d in np.nonzero(Ms[x])[0]:
                    m_=Ms[x,d]
                    if m_>Q:
                        base=U[d] if (variant=='stall' and p.get('drain_always',False)==False) else U[d]
                        start=max(U[d],Et[x]) if p.get('drain_always',False) else U[d]
                        if not p.get('drain_always',False) and U[d]<=Et[x]: continue
                        fin=start+(m_-Q)*tin
                        for c in range(4*x,4*x+4):
                            if fin>Tn[c]: stall[c,i]=fin-Tn[c]; Tn[c]=fin
        T[:,i]=Tn
    return T,stall

def score(T,lo=1,hi=100):
    dT=np.diff(T,axis=1)
    idx=np.arange(lo,hi)           # idx steps; dT[:,idx-1]... dT[:,i-1]=T[i]-T[i-1] <-> dA[:,i-1]
    a=dA[:,idx-1]; t=dT[:,idx-1]
    mape=np.mean(np.abs(t-a)/a)*100
    fin=T[:,hi-1]; pf=fin.reshape(16,4).max(1)
    af=(A[:,hi-1]-A[:,0]).reshape(16,4).max(1)
    return mape, np.mean(np.abs(pf-af)/af)*100, (pf.max()-pf.min())*1e6, (af.max()-af.min())*1e6
