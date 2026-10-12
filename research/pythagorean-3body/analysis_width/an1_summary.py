import json,glob,os,sys,math
D=sys.argv[1]
print("== final summaries")
for f in sorted(glob.glob(D+"/p3b_r3cp/*/r3c_final_*.json")):
    d=json.load(open(f))
    print(os.path.relpath(f,D), d.get("prec"), d.get("pivot"), d.get("stop_reason"), round(d.get("t_last_ok_f") or 0,3), "steps",d.get("steps_accepted"),"disc",d.get("discarded"),"wfin",d.get("wmax_final_log10"),"wall_h",round((d.get("wall_s") or 0)/3600,2))
f=D+"/p3b_r3cp/out_r4/r3c_steps_128w.jsonl"
R=[json.loads(l) for l in open(f)]
print("== 128w steps",len(R))
# width growth per step attribution
prev=None; tot=0; bins={}
for r in R:
    w=r["wmax_log10"]
    if prev is not None:
        dw=w-prev
        b=r["h_log2"]; bins.setdefault(b,[0,0.0,0]); bins[b][0]+=1; bins[b][1]+=dw; bins[b][2]+=r["n_discard"]
    prev=w
print("h_log2  steps  sum_dlog10w  discards")
for b in sorted(bins): print(b, bins[b][0], round(bins[b][1],2), bins[b][2])
# dominant component at final
import collections
dom=collections.Counter()
for r in R:
    c={"wrap":r["wR_log10"],"rem":r["wrem_log10"],"inj":r["winj_log10"]}
    dom[max(c,key=c.get)]+=1
print("dominant component counts",dom)
# trajectory of wmax by t
for T in [5,10,15,20,25,30,35,40,45,50,55,60,61]:
    rr=[r for r in R if r["t_f"]<=T]
    print("t<=",T,"wmax",rr[-1]["wmax_log10"],"steps",len(rr))
# top 15 single-step jumps
J=sorted(((R[i]["wmax_log10"]-R[i-1]["wmax_log10"],R[i]["t_f"],R[i]["h_log2"],R[i]["frame_cond"]) for i in range(1,len(R))),reverse=True)[:15]
print("top jumps (dlog10w, t, h_log2, frame_cond)");[print(j) for j in J]
