import json, csv, collections, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
g=json.load(open("results/gates_G1_G3.json")); M=json.load(open("results/metrics_all.json"))
iso=list(csv.DictReader(open("data/isolates.tsv"),delimiter="\t"))
# Fig1 isolates per year
fig,ax=plt.subplots(figsize=(7,3.2))
for sp,c in [("klebsiella","tab:red"),("ecoli","tab:blue")]:
    cnt=collections.Counter(int(r["year"]) for r in iso if r["species"]==sp); ys=sorted(y for y in cnt if y>=2008)
    ax.bar([y+(0.2 if sp=="ecoli" else -0.2) for y in ys],[cnt[y] for y in ys],0.4,color=c,label="K. pneumoniae" if sp=="klebsiella" else "E. coli")
ax.axvline(2019.5,ls="--",c="k"); ax.text(2019.6,ax.get_ylim()[1]*0.9,"date cutoff",fontsize=8)
ax.set_xlabel("collection year"); ax.set_ylabel("isolates with AST"); ax.legend(); plt.tight_layout(); plt.savefig("paper/figures/fig1_isolates_by_year.png",dpi=160); plt.close()
# Fig2 random vs date BA
cells=g["G1_cells"]; fig,axs=plt.subplots(1,3,figsize=(10,3.6),sharey=True)
for ax,m in zip(axs,["M0","M1","M2"]):
    C=[c for c in cells if c["model"]==m]
    for c in C:
        col="tab:red" if c["species"]=="klebsiella" else "tab:blue"
        ax.plot([0,1],[c["BA_random"],c["BA_date"]],"-o",c=col,alpha=.7)
    ax.set_xticks([0,1]); ax.set_xticklabels(["random split","date split"]); ax.set_title({"M0":"M0 rule","M1":"M1 logistic","M2":"M2 random forest"}[m])
axs[0].set_ylabel("balanced accuracy"); plt.tight_layout(); plt.savefig("paper/figures/fig2_random_vs_date.png",dpi=160); plt.close()
# Fig3 drops bar
fig,ax=plt.subplots(figsize=(8,3.6)); C=[c for c in cells if c["model"]=="M2"]
lab=[f'{c["species"][:4]}-{c["drug"][:8]}' for c in C]; ax.bar(range(len(C)),[c["drop"] for c in C],color=["tab:red" if c["species"]=="klebsiella" else "tab:blue" for c in C])
ax.axhline(0.02,ls="--",c="k"); ax.set_xticks(range(len(C))); ax.set_xticklabels(lab,rotation=45,ha="right",fontsize=8); ax.set_ylabel("BA drop (random - date)"); plt.tight_layout(); plt.savefig("paper/figures/fig3_drop_rf.png",dpi=160); plt.close()
# Fig4 G2
c=g["G2"]; fig,ax=plt.subplots(figsize=(5,3.2)); k=["a_unseen_determinant","b_known_determinant","c_no_determinant"]
ax.bar(["unseen\ndeterminant","known\ndeterminant","no\ndeterminant"],[c[x] for x in k],color=["tab:purple","tab:orange","tab:gray"]); ax.set_ylabel("RF false negatives (date split)"); plt.tight_layout(); plt.savefig("paper/figures/fig4_error_decomposition.png",dpi=160); plt.close()
# Fig5 G4 distances
S=list(csv.DictReader(open("results/G4_structural_alleles.tsv"),delimiter="\t"))
fig,ax=plt.subplots(figsize=(7,3.4))
for i,fam in enumerate(["blaKPC","blaSHV","blaNDM","blaIMP"]):
    for post,col in [(False,"tab:gray"),(True,"tab:red")]:
        d=[float(s["min_dist_site"]) for s in S if s["family"]==fam and s["min_dist_site"] and s["first_year"] and (int(s["first_year"])>=2020)==post]
        ax.scatter(np.full(len(d),i+(0.15 if post else -0.15))+np.random.default_rng(1).normal(0,0.04,len(d)),d,s=10,c=col,label=("first seen 2020+" if post else "first seen <2020") if i==0 else None)
ax.axhline(8,ls="--",c="k"); ax.set_xticks(range(4)); ax.set_xticklabels(["KPC (2OV5)","SHV (1SHV)","NDM (4EYL)","IMP (1DD6)"]); ax.set_ylabel("closest substitution to site (A)"); ax.legend(fontsize=8); plt.tight_layout(); plt.savefig("paper/figures/fig5_structural.png",dpi=160); plt.close()
# Fig6 date-split sensitivity by drug
fig,ax=plt.subplots(figsize=(8,3.4))
D=[m for m in M if m["split"]=="date" and m["model"]=="M2"]
ax.bar(range(len(D)),[m["sens"] for m in D],color=["tab:red" if m["species"]=="klebsiella" else "tab:blue" for m in D]); ax.set_xticks(range(len(D))); ax.set_xticklabels([f'{m["species"][:4]}-{m["drug"][:8]}' for m in D],rotation=45,ha="right",fontsize=8); ax.set_ylabel("sensitivity for R (RF, date split)"); plt.tight_layout(); plt.savefig("paper/figures/fig6_sensitivity.png",dpi=160); plt.close()
print("ok")
