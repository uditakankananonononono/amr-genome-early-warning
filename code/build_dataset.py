import csv, sys, json, re
csv.field_size_limit(10**9)
DRUGS=["ceftazidime","meropenem","ciprofloxacin","gentamicin","trimethoprim-sulfamethoxazole"]
SRC={"klebsiella":"../data/raw/PDG000000012.2532.amr.metadata.tsv","ecoli":"../data/raw/PDG000000004.6316.amr.metadata.tsv"}
out=open("data/isolates.tsv","w"); w=csv.writer(out,delimiter="\t")
w.writerow(["species","biosample","year","geo","genotypes"]+DRUGS)
n=0
for sp,p in SRC.items():
    with open(p) as f:
        r=csv.reader(f,delimiter="\t"); h=next(r); ix={k:i for i,k in enumerate(h)}
        for row in r:
            ast=row[ix["AST_phenotypes"]].strip('"')
            if ast in ("","NULL"): continue
            m=re.match(r"(\d{4})",row[ix["collection_date"]])
            if not m: continue
            ph={}
            for t in ast.split(","):
                if "=" in t:
                    k,v=t.rsplit("=",1); ph[k.strip()]=v.strip()
            if not any(d in ph for d in DRUGS): continue
            g=row[ix["AMR_genotypes"]].strip('"')
            g="" if g=="NULL" else g
            w.writerow([sp,row[ix["biosample_acc"]],m.group(1),row[ix["geo_loc_name"]].split(":")[0],g]+[ph.get(d,"") for d in DRUGS]); n+=1
print("isolates",n)
