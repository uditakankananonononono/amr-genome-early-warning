#!/usr/bin/env python3
"""amr-ew: genome-to-phenotype AMR prediction with date-aware early-warning flags.
Built from this study's findings:
 - G1: random-split accuracy overstates new-isolate accuracy by ~0.13 BA -> the tool reports the
   DATE-SPLIT accuracy of its model for each drug/species, not the random-split one.
 - G2: most new-era misses carry KNOWN determinants -> per-determinant 'miss risk' (share of
   resistant date-test isolates carrying that determinant that the model missed).
 - G3: novelty counts do not predict misses (AUROC 0.52) -> novelty is shown as information only, never as a confidence score.
 - G4: structural proximity of bla substitutions to the catalytic site is reported per allele.
Usage: amr_ew.py --species klebsiella --year 2023 --genotypes "blaKPC-2,blaSHV-11,gyrA_S83I"
"""
import argparse, csv, json, os, re, sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
H=os.path.dirname(os.path.abspath(__file__)); R=os.path.join(H,"..")
DRUGS=["ceftazidime","meropenem","ciprofloxacin","gentamicin","trimethoprim-sulfamethoxazole"]
def base(g): return re.sub(r"=.*$","",g.strip())
def load(species,cutoff):
    rows=[r for r in csv.DictReader(open(os.path.join(R,"data/isolates.tsv")),delimiter="\t") if r["species"]==species and int(r["year"])<=cutoff]
    for r in rows: r["G"]=[base(x) for x in r["genotypes"].split(",") if x]
    return rows
_CACHE={}
def main(a, emit=True):
    G=[base(x) for x in a.genotypes.split(",") if x.strip()]
    key=(a.species,a.train_until)
    if key not in _CACHE: _CACHE[key]={"rows":load(a.species,a.train_until),"models":{}}
    C=_CACHE[key]; rows=C["rows"]
    metrics=json.load(open(os.path.join(R,"results/metrics_all.json")))
    fn=list(csv.DictReader(open(os.path.join(R,"results/G2_false_negatives.tsv")),delimiter="\t"))
    struct={r["allele"]:r for r in csv.DictReader(open(os.path.join(R,"results/G4_structural_alleles.tsv")),delimiter="\t")}
    iso={r["biosample"]:r for r in csv.DictReader(open(os.path.join(R,"data/isolates.tsv")),delimiter="\t")}
    out={"species":a.species,"year":a.year,"trained_on_years_until":a.train_until,"drugs":{}}
    vocab_all=set(g for r in rows for g in r["G"])
    out["unseen_determinants"]=[g for g in G if g not in vocab_all]
    out["bla_structure"]={g:{"closest_substitution_to_site_A":struct[g]["min_dist_site"],"active_site_proximal":bool(int(struct[g]["proximal"])),"first_year_in_PDG":struct[g]["first_year"]} for g in G if g in struct}
    for d in DRUGS:
        D=[r for r in rows if r[d] in("R","S")]
        if len(D)<50: continue
        vocab=sorted({g for r in D for g in r["G"]}); vi={g:i for i,g in enumerate(vocab)}
        X=np.zeros((len(D),len(vocab)),dtype=np.int8)
        for i,r in enumerate(D):
            for g in r["G"]: X[i,vi[g]]=1
        y=np.array([int(r[d]=="R") for r in D])
        if d not in C["models"]: C["models"][d]=RandomForestClassifier(500,random_state=42,n_jobs=2).fit(X,y)
        m=C["models"][d]
        x=np.zeros((1,len(vocab)),dtype=np.int8)
        for g in G:
            if g in vi: x[0,vi[g]]=1
        p=float(m.predict_proba(x)[0,1])
        ds=[mm for mm in metrics if mm["species"]==a.species and mm["drug"]==d and mm["model"]=="M2" and mm["split"]=="date"]
        # per-determinant miss risk among resistant date-test isolates
        risk={}
        for g in G:
            Rtest=[r for r in iso.values() if r["species"]==a.species and int(r["year"])>=2020 and r[d]=="R" and g in [base(x) for x in r["genotypes"].split(",")]]
            miss=[f for f in fn if f["species"]==a.species and f["drug"]==d and g in [base(x) for x in iso[f["biosample"]]["genotypes"].split(",")]]
            if len(Rtest)>=10: risk[g]=round(len(miss)/len(Rtest),3)
        out["drugs"][d]={"call":"R" if p>=0.5 else "S","p_resistant":round(p,3),
            "honest_accuracy_on_new_isolates_BA":round(ds[0]["BA"],3) if ds else None,
            "early_warning":[f"{g}: {int(100*v)}% of resistant 2020+ isolates with this determinant were missed" for g,v in risk.items() if v>=0.3]}
    if emit: print(json.dumps(out,indent=1))
    return out
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--species",required=True,choices=["klebsiella","ecoli"])
    ap.add_argument("--year",type=int,required=True); ap.add_argument("--genotypes",required=True)
    ap.add_argument("--train-until",type=int,default=2019); main(ap.parse_args())
