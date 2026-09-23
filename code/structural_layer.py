"""G4 structural layer. Substitutions of bla alleles vs family reference mapped onto PDB structures;
distance (min heavy-atom) from substituted residue to catalytic site: Ser70 OG (class A) or Zn ions (MBL)."""
import csv, re, json, numpy as np, collections
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, PPBuilder
from scipy.stats import fisher_exact
FAM={"blaKPC":("blaKPC-2","2OV5","A","SER70"),"blaSHV":("blaSHV-1","1SHV","A","SER70"),
     "blaNDM":("blaNDM-1","4EYL","A","ZN"),"blaIMP":("blaIMP-1","1DD6","A","ZN")}
seqs={}
for rec in SeqIO.parse("../data/raw/AMRProt.fa","fasta"):
    f=rec.description.split("|")
    if len(f)>3 and f[3].startswith("bla"): seqs.setdefault(f[3],str(rec.seq).rstrip("*"))
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load("BLOSUM62"); al.open_gap_score=-10; al.extend_gap_score=-0.5; al.mode="global"
def mapping(a,b):
    A=al.align(a,b)[0]; m={}
    for (s1,e1),(s2,e2) in zip(*A.aligned):
        for k in range(e1-s1): m[s1+k]=s2+k
    return m
first={r["allele"]:int(r["first_year"]) for r in csv.DictReader(open("results/bla_allele_first_year.tsv"),delimiter="\t")}
test_novel={"blaSHV","blaNDM","blaIMP","blaKPC"}
iso=list(csv.DictReader(open("data/isolates.tsv"),delimiter="\t"))
tn=set(g for r in iso if int(r["year"])>=2020 for g in r["genotypes"].split(",") if first.get(g,0)>=2020)
fnb=set(r["biosample"] for r in csv.DictReader(open("results/G2_false_negatives.tsv"),delimiter="\t"))
fn_alleles=set(g for r in iso if r["biosample"] in fnb for g in r["genotypes"].split(",") if g in tn)
P=PDBParser(QUIET=True); out=[]
for fam,(ref,pdb,ch,site) in FAM.items():
    s=P.get_structure(pdb,f"data/structures/{pdb}.pdb"); chain=s[0][ch]
    res=[r for r in chain if r.id[0]==" " and "CA" in r]
    from Bio.PDB.Polypeptide import three_to_index,index_to_one
    pseq="".join(index_to_one(three_to_index(r.get_resname())) for r in res)
    if site=="SER70": sites=[chain[70]["OG"].coord]
    else: sites=[a.coord for a in s[0].get_atoms() if a.element=="ZN"]
    sites=np.array(sites)
    refseq=seqs[ref]; r2p=mapping(refseq,pseq)
    for allele,sq in seqs.items():
        if not re.fullmatch(fam+r"-\d+",allele) or allele==ref: continue
        m=mapping(refseq,sq); subs=[]
        for i,j in m.items():
            if refseq[i]!=sq[j]:
                d=None
                if i in r2p:
                    rr=res[r2p[i]]; d=float(min(np.linalg.norm(sites-a.coord,axis=1).min() for a in rr if a.element!="H"))
                subs.append((f"{refseq[i]}{i+1}{sq[j]}",d))
        ds=[d for _,d in subs if d is not None]
        out.append(dict(family=fam,allele=allele,first_year=first.get(allele,""),n_subs=len(subs),
            min_dist_site=min(ds) if ds else "",proximal=int(bool(ds) and min(ds)<=8.0),
            in_test_novel=int(allele in tn),in_FN=int(allele in fn_alleles),subs=";".join(f"{a}:{'' if d is None else round(d,1)}" for a,d in subs)))
with open("results/G4_structural_alleles.tsv","w") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
obs=[o for o in out if o["first_year"]!="" and o["n_subs"]>0]
post=[o for o in obs if o["first_year"]>=2020]; pre=[o for o in obs if o["first_year"]<2020]
tab=[[sum(o["proximal"] for o in post),len(post)-sum(o["proximal"] for o in post)],[sum(o["proximal"] for o in pre),len(pre)-sum(o["proximal"] for o in pre)]]
OR,p=fisher_exact(tab)
summ=dict(gate_test=dict(novel_alleles_in_test=len(tn),in_FN=sorted(fn_alleles),note="n<10 -> descriptive per GATES_LOCKED G4"),
  exploratory_post_vs_pre=dict(table_proximal_distal=tab,OR=OR,p=p,label="EXPLORATORY (not pre-registered): post-2019 vs pre-2020 alleles observed in PDG, proximal = <=8 A to site"),
  per_family={fam:dict(n=sum(o["family"]==fam for o in obs),post=sum(o["family"]==fam and o in post for o in obs)) for fam in FAM},
  test_novel_detail=[o for o in out if o["in_test_novel"]])
json.dump(summ,open("results/G4_summary.json","w"),indent=1,default=str)
print(json.dumps(summ,indent=1,default=str)[:3000])
