priority = {}
priority["Phospho[S]"] = 1e8

# pS>pT??, pS is more easy to generate modloss ions
priority["Phospho[T]"] = 1e7

priority["GG[K]"] = 1e6

# pS/pT > oM??
priority["Oxidation[M]"] = 1e5
