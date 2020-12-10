from ..spectronaut import csv as _csv

_csv._mod_dict = {
    "Carbamidomethyl[C]": "[UniMod:4]",
    "Oxidation[M]": "[UniMod:35]",
    "Phospho[S]": "[UniMod:21]",
    "Phospho[T]": "[UniMod:21]",
    "Phospho[Y]": "[UniMod:21]",
    "GlyGly[K]": "[UniMod:121]",
}

DIANN_CSV = _csv.SPN_CSV