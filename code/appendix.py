import json,csv,collections,re
M=json.load(open("results/metrics_all.json")); fn=list(csv.DictReader(open("results/G2_false_negatives.tsv"),delimiter="\t"))
iso=list(csv.DictReader(open("data/isolates.tsv"),delimiter="\t"))
L=['<div style="page-break-after: always;"></div>',"","# Appendix A: all 60 model-split cells","",
"| Species | Drug | Split | Model | n test | n R | BA (95% CI) | MCC | Sens. |","|---|---|---|---|---|---|---|---|---|"]
for m in M: L.append(f'| {m["species"]} | {m["drug"]} | {m["split"]} | {m["model"]} | {m["n_test"]} | {m["n_R"]} | {m["BA"]:.3f} ({m["BA_CI"][0]:.3f}-{m["BA_CI"][1]:.3f}) | {m["MCC"]:.3f} | {m["sens"]:.3f} |')
L+=["","# Appendix B: error anatomy by drug and species","","| Species | Drug | (a) unseen | (b) known | (c) none | total |","|---|---|---|---|---|---|"]
c=collections.defaultdict(collections.Counter)
for f in fn: c[(f["species"],f["drug"])][f["category"]]+=1
for k,v in sorted(c.items()): L.append(f'| {k[0]} | {k[1]} | {v["a_unseen_determinant"]} | {v["b_known_determinant"]} | {v["c_no_determinant"]} | {sum(v.values())} |')
L+=["","# Appendix C: known determinants most often present in missed resistant isolates","","Share = missed resistant 2020+ carriers / all resistant 2020+ carriers (random forest, date split). Determinants with at least 20 resistant carriers; top 30 by count of misses.","",
"| Species | Drug | Determinant | Missed | Resistant carriers | Share missed |","|---|---|---|---|---|---|"]
base=lambda g: re.sub(r"=.*$","",g)
G={r["biosample"]:set(base(x) for x in r["genotypes"].split(",") if x) for r in iso}
rows=[]
for sp in ["klebsiella","ecoli"]:
  for d in ["ceftazidime","meropenem","ciprofloxacin","gentamicin","trimethoprim-sulfamethoxazole"]:
    R=[r for r in iso if r["species"]==sp and int(r["year"])>=2020 and r[d]=="R"]
    mis=collections.Counter(g for f in fn if f["species"]==sp and f["drug"]==d for g in G[f["biosample"]])
    car=collections.Counter(g for r in R for g in G[r["biosample"]])
    for g,n in car.items():
      if n>=20: rows.append((mis[g],sp,d,g,n))
for m,sp,d,g,n in sorted(rows,reverse=True)[:30]: L.append(f"| {sp} | {d} | {g} | {m} | {n} | {m/n:.2f} |")
L+=["","Many top entries are intrinsic or near-universal genes (e.g. fosA, oqxA/B, blaSHV variants in K. pneumoniae); their presence in missed isolates is expected and not itself causal. Acquired determinants in this table are the candidates for determinant-level drift monitoring.","",
"# Appendix D: gate ledger","","| Gate | Criterion (locked) | Outcome |","|---|---|---|",
"| G0 | >=3,000 isolates; each test cell >=100 with >=20 R | PASS (12,841; all 10 cells) |",
"| G1 | mean ML BA drop >=0.02, CI excludes 0 | PASS (0.134, CI 0.109-0.158) |",
"| G2 | decompose all date-split FNs | DONE (41 / 869 / 371) |",
"| G3 | novelty AUROC >=0.65, CI lower >0.5 | FAIL (0.522) - kept as negative |",
"| G4 | Fisher test if n>=10, else descriptive | DESCRIPTIVE (n=5); exploratory null |",
"| G5 | CLI runs on held-out date-split isolates | PASS (77 calls, accuracy 0.870) |"]
open("paper/appendix.md","w").write("\n".join(L)+"\n")
