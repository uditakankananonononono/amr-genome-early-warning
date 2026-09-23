"""G5: run the tool's model path on held-out date-split isolates (trained <=2019) and check it matches locked M2 metrics; also measure whether early-warning flags concentrate misses."""
import csv, json, re, sys, numpy as np, subprocess, os, random
sys.path.insert(0,os.path.dirname(__file__))
random.seed(42)
import amr_ew, argparse
iso=[r for r in csv.DictReader(open("data/isolates.tsv"),delimiter="\t") if int(r["year"])>=2020]
sample=random.sample([r for r in iso if r["species"]=="klebsiella"],15)+random.sample([r for r in iso if r["species"]=="ecoli"],15)
ok=0; rows=[]
for r in sample:
    o=amr_ew.main(argparse.Namespace(species=r["species"],year=int(r["year"]),genotypes=r["genotypes"] or "none",train_until=2019),emit=False)
    for d,v in o["drugs"].items():
        if r[d] in("R","S"): rows.append((r["species"],d,r[d],v["call"],int(bool(v["early_warning"]))))
acc=np.mean([t==c for _,_,t,c,_ in rows])
missed=[x for x in rows if x[2]=="R" and x[3]=="S"]; caught=[x for x in rows if x[2]=="R" and x[3]=="R"]
res=dict(n_isolates=len(sample),n_calls=len(rows),accuracy=float(acc),
    flag_rate_in_missed_R=float(np.mean([x[4] for x in missed])) if missed else None,n_missed=len(missed),
    flag_rate_in_caught_R=float(np.mean([x[4] for x in caught])) if caught else None,n_caught=len(caught),
    note="30 random 2020+ isolates (seed 42), tool trained on <=2019. Flags use per-determinant miss rates computed on the SAME date-test set, so the flag rates here are in-sample and optimistic; reported for transparency.")
json.dump(res,open("results/G5_tool_validation.json","w"),indent=1); print(json.dumps(res,indent=1))
