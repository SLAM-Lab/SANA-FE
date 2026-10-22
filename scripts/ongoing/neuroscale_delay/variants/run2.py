from sim2 import *
rep('V2 no stall',sim2(P,True,True)[0])
for H in (0,1,2,7):
    T,s=sim2(dict(P,Q=50,Hmax=H),True,True,'hard'); rep('V3 hist=%d Q=50'%H,T,s)
# distribution of stall magnitude actual vs model
T2,_=sim2(P,True,True); T3,s3=sim2(dict(P,Q=50),True,True,'hard')
# per-tile finish at several steps
for st in (25,50,75,99):
    a=(A[:,st]-A[:,0]).reshape(16,4).max(1); t2=T2[:,st].reshape(16,4).max(1); t3=T3[:,st].reshape(16,4).max(1)
    print('step',st,'mean|err| V2 %.2f%%  V3 %.2f%%  (us: actsim span %.0f)'%(100*np.mean(abs(t2-a)/a),100*np.mean(abs(t3-a)/a),(a.max()-a.min())*1e6))
# windowed MAPE (10-step) per core
def wm(T,w=10):
    a=A-A[:,[0]]; e=[]
    for s in range(w,100,w): e.append(np.abs((T[:,s]-T[:,s-w])-(a[:,s]-a[:,s-w]))/(a[:,s]-a[:,s-w]))
    return 100*np.mean(e)
print('10-step window error V2 %.2f%% V3 %.2f%%'%(wm(T2),wm(T3)))
# aggregate stall contribution
print('mean stall per core-step (us):',s3.mean()*1e6,' vs mean step',np.diff(T3,axis=1).mean()*1e6)
