from harness import *
nax=np.load(R+'nax_static.npy'); cum=np.cumsum(mi,axis=1)>0
sync=n_fi+n_fo
P=dict(F=14.0e-6+50*7.6e-9,t_scan=3.8e-9,t_ev=4.7e-9,t_in=11.2e-9,t_out=3.0e-9,t_pkt=3.17e-9,W=14.0e-6,Q=50,tdrain=11.2e-9,Hmax=7,p_exp=1.0)
Ms_all=M.reshape(NS,16,4,64).sum(2)    # [i, srctile, dstcore]
def sim2(p,scan=False,port=False,stall=None):
    T=np.zeros((NC,NS)); Sall=np.zeros((NC,NS)); Wend=np.zeros((NC,NS)); st=np.zeros((NC,NS))
    for i in range(1,NS):
        prev=T[:,i-1]; S=prev.copy()
        for c in range(NC):
            if len(FI[c]): S[c]=max(S[c],prev[FI[c]].max())
            if len(FO[c]) and i-TADV>=0: S[c]=max(S[c],T[FO[c],i-TADV].max())
        acc=p['t_ev']*ev[:,i-1]
        L=p['F']+acc+p['t_in']*mi[:,i]
        if scan: L=L+p['t_scan']*nax*cum[:,i-1]
        L=L+(p['t_pkt']*tilesum(mo[:,i]+sync) if port else p['t_out']*mo[:,i])
        Tn=S+L
        Sall[:,i]=S; Wend[:,i]=S+acc+p['W']
        if stall:
            E=Wend[:,i].reshape(16,4).max(1)
            out=np.zeros(16)
            for x in range(16):
                for d in np.nonzero(Ms_all[i,x])[0]:
                    m_=Ms_all[i,x,d]
                    if m_<=p['Q']: continue
                    for j in range(max(1,i-p['Hmax']),i+1):
                        a=Sall[d,j]; b=Wend[d,j]
                        if j==i and d: pass
                        if a<=E[x]<b:
                            out[x]=max(out[x],b+(m_-p['Q'])*p['tdrain']); break
            for x in range(16):
                if out[x]>0:
                    for c in range(4*x,4*x+4):
                        if out[x]>Tn[c]:
                            add=out[x]-Tn[c]
                            if stall=='expected': add*=p['p_exp']
                            st[c,i]=add; Tn[c]+=add
        T[:,i]=Tn
    return T,st
def rep(name,T,st=None):
    a=score(T,1,100); b=score(T,50,100)
    print('%-44s MAPE %5.2f%% (test %5.2f%%) pop-finish err %5.2f%% spread %4.0fus (actsim %.0f) stalls %s'%(name,a[0],b[0],a[1],a[2],a[3],'-' if st is None else int((st>0).sum())))
if __name__=='__main__':
    rep('V0 YAML additive',sim2(P)[0])
    rep('V1 +scan (num_ax*3.8ns)',sim2(P,scan=True)[0])
    rep('V2 +tile port, sync packets',sim2(P,True,True)[0])
    for Q in (30,40,50,66):
        T,s=sim2(dict(P,Q=Q),True,True,'hard'); rep('V3 +dest-window stall Q=%d'%Q,T,s)
    for pe in (0.35,0.5):
        T,s=sim2(dict(P,Q=40,p_exp=pe),True,True,'expected'); rep('V3e expected-stall Q=40 p=%.2f'%pe,T,s)
