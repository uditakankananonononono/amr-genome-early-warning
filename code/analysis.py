"""G1-G3 per GATES_LOCKED.md. Seed 42."""
import csv, json, re, numpy as np
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, recall_score, roc_auc_score
from sklearn.model_selection import train_test_split
rng=np.random.default_rng(42)
DRUGS=["ceftazidime","meropenem","ciprofloxacin","gentamicin","trimethoprim-sulfamethoxazole"]
# drug -> determinant subclass tokens (M0 rule). SXT needs trimethoprim AND sulfonamide determinant.
RULE={"ceftazidime":[{"CEPHALOSPORIN","CARBAPENEM"}],"meropenem":[{"CARBAPENEM"}],
      "ciprofloxacin":[{"QUINOLONE"}],"gentamicin":[{"GENTAMICIN"}],
      "trimethoprim-sulfamethoxazole":[{"TRIMETHOPRIM"},{"SULFONAMIDE"}]}
sub={}
for r in csv.DictReader(open("data/ReferenceGeneCatalog.txt"),delimiter="\t"):
    s=set(r["subclass"].split("/"))|({r["class"]} if r["class"] else set())
    for k in (r["allele"],r["gene_family"]):
        if k: sub.setdefault(k,set()).update(s)
for f,col in [("data/AMRProt-mutation.tsv","standard_mutation_symbol"),("data/AMR_DNA-Escherichia.tsv","standard_mutation_symbol"),("data/AMR_DNA-Klebsiella_pneumoniae.tsv","standard_mutation_symbol")]:
    rd=csv.reader(open(f),delimiter="\t"); h=next(rd); h[0]=h[0].lstrip("#"); ix={k:i for i,k in enumerate(h)}
    for r in rd: sub.setdefault(r[ix[col]],set()).update(r[ix["subclass"]].split("/"))
def base(g): return re.sub(r"=.*$","",g)
def subclasses(g):
    b=base(g)
    if b in sub: return sub[b]
    fam=re.sub(r"-\d+.*$","",b)
    return sub.get(fam,set())
rows=list(csv.DictReader(open("data/isolates.tsv"),delimiter="\t"))
for r in rows: r["G"]=[base(x) for x in r["genotypes"].split(",") if x]; r["Y"]=int(r["year"])
def m0(r,d): return int(all(any(subclasses(g)&need for g in r["G"]) for need in RULE[d]))
def boot(yt,a,b,n=1000):
    yt=np.array(yt);a=np.array(a);b=np.array(b);o=[]
    for _ in range(n):
        i=rng.integers(0,len(yt),len(yt))
        if len(set(yt[i]))<2: continue
        o.append(balanced_accuracy_score(yt[i],a[i])-balanced_accuracy_score(yt[i],b[i]))
    return np.percentile(o,[2.5,97.5]).tolist()
def bootci(yt,p,n=1000):
    yt=np.array(yt);p=np.array(p);o=[]
    for _ in range(n):
        i=rng.integers(0,len(yt),len(yt))
        if len(set(yt[i]))<2: continue
        o.append(balanced_accuracy_score(yt[i],p[i]))
    return np.percentile(o,[2.5,97.5]).tolist()
res=[];fn_rows=[];ew=[]
for sp in ["klebsiella","ecoli"]:
  for d in DRUGS:
    D=[r for r in rows if r["species"]==sp and r[d] in("R","S")]
    y=np.array([int(r[d]=="R") for r in D])
    tr_d=[i for i,r in enumerate(D) if r["Y"]<=2019]; te_d=[i for i,r in enumerate(D) if r["Y"]>=2020]
    tr_r,te_r=train_test_split(np.arange(len(D)),test_size=0.2,stratify=y,random_state=42)
    for split,tr,te in [("random",tr_r,te_r),("date",tr_d,te_d)]:
        vocab=sorted({g for i in tr for g in D[i]["G"]}); vi={g:k for k,g in enumerate(vocab)}
        def X(idx):
            M=np.zeros((len(idx),len(vocab)),dtype=np.int8)
            for a,i in enumerate(idx):
                for g in D[i]["G"]:
                    if g in vi: M[a,vi[g]]=1
            return M
        Xtr,Xte=X(tr),X(te); ytr,yte=y[tr],y[te]
        preds={"M0":np.array([m0(D[i],d) for i in te]),
               "M1":LogisticRegression(C=1.0,max_iter=2000).fit(Xtr,ytr).predict(Xte),
               "M2":RandomForestClassifier(500,random_state=42,n_jobs=2).fit(Xtr,ytr).predict(Xte)}
        for m,p in preds.items():
            res.append(dict(species=sp,drug=d,split=split,model=m,n_test=len(te),n_R=int(yte.sum()),
                BA=balanced_accuracy_score(yte,p),BA_CI=bootci(yte,p),MCC=matthews_corrcoef(yte,p),sens=recall_score(yte,p),_y=yte.tolist(),_p=p.tolist()))
        if split=="date":
            trset=set(vocab)
            p2=preds["M2"]
            for a,i in enumerate(te):
                r=D[i]; unseen=[g for g in r["G"] if g not in trset]
                unseen_rel=[g for g in unseen if any(subclasses(g)&need for need in RULE[d])]
                known_rel=[g for g in r["G"] if g in trset and any(subclasses(g)&need for need in RULE[d])]
                if yte[a]==1:
                    ew.append(dict(species=sp,drug=d,fn=int(p2[a]==0),novelty=len(unseen),novelty_rel=len(unseen_rel)))
                if yte[a]==1 and p2[a]==0:
                    cat="a_unseen_determinant" if unseen_rel else ("b_known_determinant" if known_rel else "c_no_determinant")
                    fn_rows.append(dict(species=sp,drug=d,biosample=r["biosample"],year=r["Y"],category=cat,unseen_relevant=";".join(unseen_rel),known_relevant=";".join(known_rel)))
# G1
g1=[]
for sp in ["klebsiella","ecoli"]:
  for d in DRUGS:
    for m in ["M0","M1","M2"]:
      R=[x for x in res if x["species"]==sp and x["drug"]==d and x["model"]==m]
      rr=[x for x in R if x["split"]=="random"][0]; dd=[x for x in R if x["split"]=="date"][0]
      g1.append(dict(species=sp,drug=d,model=m,BA_random=rr["BA"],BA_date=dd["BA"],drop=rr["BA"]-dd["BA"]))
drops=np.array([x["drop"] for x in g1 if x["model"]!="M0"])
bs=[rng.choice(drops,len(drops)).mean() for _ in range(1000)]
g1sum=dict(mean_drop_ML=float(drops.mean()),CI=np.percentile(bs,[2.5,97.5]).tolist(),n_cells=len(drops),
           PASS=bool(drops.mean()>=0.02 and np.percentile(bs,2.5)>0))
# G3
yfn=[e["fn"] for e in ew]; sc=[e["novelty"]+2*e["novelty_rel"] for e in ew]
auc=roc_auc_score(yfn,sc); ab=[]
yfn_a=np.array(yfn); sc_a=np.array(sc)
for _ in range(1000):
    i=rng.integers(0,len(yfn_a),len(yfn_a))
    if len(set(yfn_a[i]))<2: continue
    ab.append(roc_auc_score(yfn_a[i],sc_a[i]))
g3=dict(AUROC=auc,CI=np.percentile(ab,[2.5,97.5]).tolist(),n_R=len(ew),n_FN=int(sum(yfn)),PASS=bool(auc>=0.65 and np.percentile(ab,2.5)>0.5),
        score="novelty = #determinants unseen before cutoff + 2 x #unseen determinants of the drug's class")
cats=defaultdict(int)
for f in fn_rows: cats[f["category"]]+=1
json.dump([{k:v for k,v in x.items() if not k.startswith("_")} for x in res],open("results/metrics_all.json","w"),indent=1)
json.dump(dict(G1=g1sum,G1_cells=g1,G2=dict(cats),G3=g3),open("results/gates_G1_G3.json","w"),indent=1)
with open("results/G2_false_negatives.tsv","w") as f:
    w=csv.DictWriter(f,fieldnames=list(fn_rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(fn_rows)
print(json.dumps(g1sum)); print(dict(cats)); print(json.dumps(g3))
for x in g1: print(f'{x["species"]:10s} {x["drug"][:14]:14s} {x["model"]} rand {x["BA_random"]:.3f} date {x["BA_date"]:.3f} drop {x["drop"]:+.3f}')
