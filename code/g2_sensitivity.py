"""POST-HOC sensitivity for G2 (not pre-registered): re-classify RF date-split false negatives after
removing near-universal (likely intrinsic) determinants, defined as present in >=80% of that species' isolates."""
import csv,re,json,collections
exec(open("code/analysis.py").read().split("rows=list(")[0])
iso=list(csv.DictReader(open("data/isolates.tsv"),delimiter="\t"))
fn=list(csv.DictReader(open("results/G2_false_negatives.tsv"),delimiter="\t"))
G={r["biosample"]:[base(x) for x in r["genotypes"].split(",") if x] for r in iso}
univ={}
for sp in ["klebsiella","ecoli"]:
    S=[r for r in iso if r["species"]==sp]; c=collections.Counter(g for r in S for g in set(G[r["biosample"]]))
    univ[sp]={g for g,n in c.items() if n/len(S)>=0.8}
cats=collections.Counter()
for f in fn:
    unseen=set(f["unseen_relevant"].split(";"))-{""}
    known=set(f["known_relevant"].split(";"))-{""}-univ[f["species"]]
    cats["a_unseen_determinant" if unseen else ("b_known_acquired_determinant" if known else "c_no_acquired_determinant")]+=1
out=dict(label="POST-HOC sensitivity, not pre-registered",intrinsic_proxy={k:sorted(v) for k,v in univ.items()},categories=dict(cats),n=len(fn))
json.dump(out,open("results/G2_sensitivity_intrinsic.json","w"),indent=1); print(json.dumps(out))
