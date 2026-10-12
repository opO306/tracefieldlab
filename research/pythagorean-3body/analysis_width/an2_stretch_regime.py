import json,gzip,sys,math
import numpy as np
D=sys.argv[1]
R=[json.loads(l) for l in open(D+"/p3b_r3cp/out_r4/r3c_steps_128w.jsonl")]
ref=json.load(open(D+"/p3b_an/results/ref_B_grid.json"))
st=np.array([[float(x) for x in s[:6]] for s in ref["states"]]); dt=0.01
def rmin(t):
    k=min(int(round(t/dt)),len(st)-1); p=st[k].reshape(3,2)
    return min(np.linalg.norm(p[0]-p[1]),np.linalg.norm(p[0]-p[2]),np.linalg.norm(p[1]-p[2]))
# true stretching via QR accumulation of step Jacobians
Q=np.eye(12); lam=np.zeros(12); out=[]
with gzip.open(D+"/p3b_r3cp/out_r4/r3c_jac_128w.jsonl.gz","rt") as f:
    for line in f:
        j=json.loads(line); M=np.array(j["M"]).reshape(12,12)
        Z=M@Q; Q,Rr=np.linalg.qr(Z); lam+=np.log10(np.abs(np.diag(Rr)))
        out.append((j["i"],j["t_f"],lam[0]))
print("jac records",len(out),"steps",len(R))
s=dict((i,(t,l)) for i,t,l in out)
w0=R[0]["wmax_log10"]
print(" t      log10w-w0  log10stretch  excess   rmin")
for T in [1,5,10,15,20,25,30,35,40,45,50,55,58,60,61]:
    r=[x for x in R if x["t_f"]<=T][-1]; t,l=s[r["i"]]
    print(f"{T:4} {r['wmax_log10']-w0:10.2f} {l:12.2f} {r['wmax_log10']-w0-l:8.2f}  {rmin(t):.4f}")
# growth per unit time by regime of rmin
bins={}
for a,b in zip(R,R[1:]):
    rm=rmin(b["t_f"]); key="r<0.05" if rm<0.05 else ("r<0.2" if rm<0.2 else ("r<0.5" if rm<0.5 else "r>=0.5"))
    dw=b["wmax_log10"]-a["wmax_log10"]; dl=s[b["i"]][1]-s[a["i"]][1]; h=b["t_f"]-a["t_f"]
    e=bins.setdefault(key,[0,0.0,0.0,0.0,0]); e[0]+=1;e[1]+=h;e[2]+=dw;e[3]+=dl;e[4]+=b["n_discard"]
print("regime steps time dlog10w dlog10stretch discards  w/time")
for k,e in bins.items(): print(k,e[0],round(e[1],2),round(e[2],2),round(e[3],2),e[4],round(e[2]/e[1],3))
