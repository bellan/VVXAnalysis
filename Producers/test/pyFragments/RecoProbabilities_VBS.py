
AJetsProdDecProbabilities_MCFM_JECNominal = [
   "Name:JJEW_BKG_MCFM_JECNominal Process:bkgZZ Production:JJEW MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJVBF_BKG_MCFM_JECNominal Process:bkgZZ Production:JJVBF MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
   "Name:JJQCD_BKG_MCFM_JECNominal Process:bkgZZ Production:JJQCD MatrixElement:MCFM Cluster:J2JECNominal DefaultME:-1",
]
AJetsProdDecProbabilities_MCFM_JECUp = [theME.replace("JECNominal", "JECUp") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]
AJetsProdDecProbabilities_MCFM_JECDn = [theME.replace("JECNominal", "JECDn") for theME in AJetsProdDecProbabilities_MCFM_JECNominal]



# Construct the final list
theRecoProbabilities = []

theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECNominal)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECUp)
theRecoProbabilities.extend(AJetsProdDecProbabilities_MCFM_JECDn)

# Append final list
for name in (
             "ZZCand",          "ZZTree",
             "ZLLCand",         "CRZLLTree",
             "ZZCandlooseEle",  "ZZTreelooseEle",
             "ZLLCandlooseEle", "CRZLLTreelooseEle",
             "ZLLCandZ1RSE",    "CRZLLTreeZ1RSE",
             "ZZCandtle",       "ZZTreetle",
             "ZLLCandtle",      "CRZLLTreetle",
            ):
    if hasattr(process, name):
        getattr(process, name).recoProbabilities.extend(theRecoProbabilities)
