"""Example molecules: global-health drugs, as (SMILES, name).

FAMILIES share a Bemis-Murcko scaffold, so "Group by scaffold" has analogues to align.
SINGLES are structurally unrelated. Checked against their molecular formulas.
"""

import random

FAMILIES = [
    [  # 4-aminoquinolines
        ("CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "Chloroquine"),
        ("CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "Hydroxychloroquine"),
        ("CCNCCCC(C)Nc1ccnc2cc(Cl)ccc12", "Desethylchloroquine"),
    ],
    [  # artemisinin derivatives (lactol ethers and esters)
        ("C[C@@H]1CC[C@H]2[C@H]([C@H](O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)O)C", "Dihydroartemisinin"),
        ("C[C@@H]1CC[C@H]2[C@H]([C@H](O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)OC)C", "Artemether"),
        ("C[C@@H]1CC[C@H]2[C@H]([C@H](O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)OC(=O)CCC(=O)O)C", "Artesunate"),
    ],
    [  # penicillins
        ("CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", "Penicillin G"),
        ("CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O", "Ampicillin"),
        ("CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", "Amoxicillin"),
    ],
    [  # fluoroquinolones
        ("O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", "Ciprofloxacin"),
        ("CCN1CCN(CC1)c1cc2n(C3CC3)cc(C(=O)O)c(=O)c2cc1F", "Enrofloxacin"),
    ],
    [  # 5-nitroimidazoles
        ("Cc1ncc([N+](=O)[O-])n1CCO", "Metronidazole"),
        ("CCS(=O)(=O)CCn1c([N+](=O)[O-])cnc1C", "Tinidazole"),
        ("CC(O)Cn1c([N+](=O)[O-])cnc1C", "Secnidazole"),
        ("Cc1ncc([N+](=O)[O-])n1CC(O)CCl", "Ornidazole"),
    ],
    [  # pyridine antitubercular agents
        ("NNC(=O)c1ccncc1", "Isoniazid"),
        ("CCc1cc(C(N)=S)ccn1", "Ethionamide"),
    ],
]

SINGLES = [
    ("C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C", "Artemisinin"),
    ("O[C@H](c1cc(C(F)(F)F)nc2c(C(F)(F)F)cccc12)[C@@H]1CCCCN1 |&1:1,20|", "rac-Mefloquine"),
    ("CCc1nc(N)nc(N)c1-c1ccc(Cl)cc1", "Pyrimethamine"),
    ("COc1cc(NC(C)CCCN)c2ncccc2c1", "Primaquine"),
    ("COc1cc(NC(C)CCCN)c2nc(OC)cc(C)c2c1Oc1cccc(C(F)(F)F)c1", "Tafenoquine"),
    ("COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", "Trimethoprim"),
    ("CCCCN(CCCC)CC(O)c1cc(Cl)cc2/C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12", "Lumefantrine"),
    ("O=C1C(O)=C([C@H]2CC[C@H](c3ccc(Cl)cc3)CC2)C(=O)c2ccccc12", "Atovaquone"),
    ("O=C(Cn1ccnc1[N+](=O)[O-])NCc1ccccc1", "Benznidazole"),
    ("CC1CS(=O)(=O)CCN1/N=C/c1ccc(o1)[N+](=O)[O-]", "Nifurtimox"),
    ("Cc1ncc([N+](=O)[O-])n1COc1ccc(SC)cc1", "Fexinidazole"),
    ("CCCCCCCCCCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C", "Miltefosine"),
    ("O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1", "Praziquantel"),
    ("NC(=O)c1cnccn1", "Pyrazinamide"),
    ("CC[C@@H](CO)NCCN[C@@H](CC)CO", "Ethambutol"),
    ("CC(=O)NC[C@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1", "Linezolid"),
    ("COc1nc2ccc(Br)cc2cc1[C@@H](c1ccccc1)[C@@](O)(CCN(C)C)c1cccc2ccccc12", "Bedaquiline"),
    ("[O-][N+](=O)c1cn2C[C@H](OCc3ccc(OC(F)(F)F)cc3)COc2n1", "Pretomanid"),
    ("Nc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1", "Dapsone"),
    ("Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "Sulfamethoxazole"),
    ("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "Fluconazole"),
    ("COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1", "Quinine"),
    ("CCN(CC)Cc1cc(Nc2ccnc3cc(Cl)ccc23)ccc1O", "Amodiaquine"),
    ("Clc1ccc2c(N3CCN(CCCN4CCN(c5ccnc6cc(Cl)ccc56)CC4)CC3)ccnc2c1", "Piperaquine"),
    ("CC(C)N=C(N)N=C(N)Nc1ccc(Cl)cc1", "Proguanil"),
    ("COc1ncnc(NS(=O)(=O)c2ccc(N)cc2)c1OC", "Sulfadoxine"),
    ("Nc1ccc(S(=O)(=O)Nc2ncccn2)cc1", "Sulfadiazine"),
    ("CC(C)N=c1cc2n(-c3ccc(Cl)cc3)c3ccccc3nc-2cc1Nc1ccc(Cl)cc1", "Clofazimine"),
    ("CCCSc1ccc2[nH]c(NC(=O)OC)nc2c1", "Albendazole"),
    ("COC(=O)Nc1nc2cc(C(=O)c3ccccc3)ccc2[nH]1", "Mebendazole"),
    ("O=C(Nc1ccc([N+](=O)[O-])cc1Cl)c1cc(Cl)ccc1O", "Niclosamide"),
    ("CCn1cc(C(=O)O)c(=O)c2cc(F)c(N3CCNCC3)cc21", "Norfloxacin"),
    ("COc1c(N2C[C@@H]3CCCN[C@@H]3C2)c(F)cc2c(=O)c(C(=O)O)cn(C3CC3)c12", "Moxifloxacin"),
    ("Cc1onc(-c2ccccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12", "Cloxacillin"),
    ("N=C(N)c1ccc(OCCCCCOc2ccc(C(N)=N)cc2)cc1", "Pentamidine"),
    ("NCCCC(N)(C(F)F)C(=O)O", "Eflornithine"),
    ("Nc1ccc(C(=O)O)c(O)c1", "Aminosalicylic acid"),
    ("N[C@@H]1CONC1=O", "Cycloserine"),
]

POOL = [m for fam in FAMILIES for m in fam] + SINGLES
SIZES = (20, 50)  # random examples hold this many molecules


def as_text(pairs):
    return "".join(f"{smi} {name}\n" for smi, name in pairs)


def shuffled(n=None, rng=random):
    """Analogue families plus unrelated singles, in random order.

    ``n`` defaults to a random size between 20 and 50 molecules. Three families are always
    included, so grouping by scaffold has analogues to align.
    """
    n = rng.randint(*SIZES) if n is None else n
    if not 1 <= n <= len(POOL):
        raise ValueError(f"n must be between 1 and {len(POOL)}, not {n}")
    families = rng.sample(FAMILIES, len(FAMILIES))
    picked, k = [], 0
    while k < len(families) and (k < 3 or n - len(picked) > len(SINGLES)):
        picked += families[k]
        k += 1
    picked = picked[:n] + rng.sample(SINGLES, max(0, n - len(picked)))
    rng.shuffle(picked)
    return picked


DEFAULT = shuffled(24, random.Random(7))  # the example shown first: fixed, 24 molecules
