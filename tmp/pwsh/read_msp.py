import sys
input = sys.argv[1]
output = sys.argv[2]
mod_dict = {
    "Acetyl": "Acetyl[ProteinN-term]",
    "Phospho": "Phospho",
    "Oxidation": "Oxidation",
    "CAM": "Carbamidomethyl",
}
with open(input) as f:
    modseq = set()
    while True:
        line = f.readline()
        if not line: break
        if not line.startswith("Comment"): continue
        items = line.split(" ")
        modlist = []
        mod_in_dict = True
        for item in items:
            if item.startswith("Mods"):
                mods = item.split("(")
                for mod in mods[1:]:
                    mod = mod[:-1].split(",")
                    if mod[-1] == "Acetyl":
                        modlist.append((0, "Acetyl[ProteinN-term]"))
                    elif mod[-1] in mod_dict:
                        modlist.append((int(mod[0])+1, "%s[%s]"%(mod_dict[mod[2]], mod[1])))
                    else:
                        mod_in_dict = False
                        break
                if not mod_in_dict: break
            elif item.startswith("Fullname="):
                seq = item.split('.')[1]
            elif item.startswith("Charge="):
                charge = item[item.find("=")+1:]
            elif item.startswith("Protein="):
                protein = item[item.find("=")+2:]
                idx = protein.rfind('"')
                if idx != -1:
                    protein = protein[:idx]
        if not mod_in_dict: continue
        # idx = seq.find("C")
        # while idx != -1:
            # modlist.append((idx+1, "Carbamidomethyl[C]"))
            # idx = seq.find("C",idx+1)
        modlist.sort(key = lambda x: x[0])
        mod = ";".join(["%d,%s"%(site, m) for site, m in modlist])
        modseq.add((seq, mod, charge, protein))
with open(output,"w") as f:
    f.write("peptide\tmodinfo\tcharge\tprotein\n")
    for item in modseq: f.write("\t".join(item)+"\n")
        