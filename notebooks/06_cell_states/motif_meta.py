clustername_mapping = {
    "Average_1139":"POU-like",
    "Average_1034":"P53-like",
    "Average_719":"SIX/ZNF32",
    "Average_1111":"TEAD-like",
    "Average_1147":"ETS-like",
    "Average_1069":"RFX", 
    "Average_1121":"BACH/JUN-like",
    "Average_1021":"CEPB-like",
    "Average_1047":"NFI-monomer",
    "Average_459":"NFI-dimer",
    "Average_1109":"SOX-like",
    "Average_608":"CTCF", 
    "Average_1142":"ATF/CREB/RXR", 
    "Average_1146":"KLF/SP",
    "Average_1135":"FOX-like",
    "Average_1149":"GRHL-like",
    "Average_1154":"E-box (CACCTG)",
    'Average_1113':'PAX/CUX/ONEC',
    'Average_245':'HNF4',
    'Average_1090':'GATA/TAL/SIX',
    'Average_1099':'STAT/IRF',
    'Average_1038':'TYY/TAF1',
    'Average_1018':'MAF-like',
    "Average_1080":'E-box (CACGTG)',
    "Average_1085":"RUNX",
    "Average_870":"MEF2",
    "Average_563":"Steroid",
    "REST.H12CORE.0.P.B":"REST",
    "Average_976":"CDX",
    'Average_1151':'KMT2A/PATZ1',
    'Average_1138': 'HNF1/SATB1',
    'ZN362.H12CORE.0.P.C':'ZN362'
}

bad_motifs = set(['Average_1151', # low-info long GC noise
                "ZBT40.H12CORE.0.P.B",# "matches" look very much like PolyA signal
                "Average_1138", # low-info AT-rich motif
                "ZN132.H12CORE.0.P.C", # low-info long
                "ZN362.H12CORE.0.P.C", # just a bunch of A's
                'Average_1117', # ZincFinger
                'ZFP14.H12CORE.0.P.C', # ZincFinger
                'Average_1145', # ZingFincer
                'ZN787.H12CORE.0.M.C', # ZingFinger
                'ZN821.H12CORE.0.SM.B', # ZincFinger
                'ZN540.H12CORE.0.P.C', # ZincFinger
                ])


