Samples to be used for the analysis
-----------------------------------------------
-----------------------------------------------

All samples processed in the analysis work-flow are listed here VVXAnalysis/Producers/python/samples_8TeV.csv, however, not all samples are mutually exclusive.

We can devide the them into two possible sets that contains (almost) totally mutually exclusive samples.


The MadGraph set
-----------------------------------------------

- qq,qg,gg -> 4l +0,1,2 j (MadGraph). LO for the zero jet bin, NLO for the others
   - ZZJetsTo4L 
- gg -> 4l (MCFM). LO. Nickname: gg_4l.
   - ggTo4e_SMHContinInterf-MCFM67_H125.6
   - ggTo4mu_SMHContinInterf-MCFM67_H125.6
   - ggTo2e2mu_SMHContinInterf-MCFM67_H125.6
- qq -> 4l + 2 j (Phantom samples). LO. Nickname: qq_4l2j.
   - ZZTo4eJJ_SMHContinInterf_H125.6
   - ZZTo2e2muJJ_SMHContinInterf_H125.6
   - ZZTo4muJJ_SMHContinInterf_H125.6
- Higgs samples
   - powheg15jhuGenV3H126
   - VBFH126
   - ttH126
   - ZH
   - WH
- Dibosons
   - WZ
   - WWJets
   - WGToLNuG
- Tribosons
   - WZZJets (or WZZ_aMCatNLO)
   - WWWJets
   - WWZJets
   - ZZZJets
   - WWGJets
- ttbar
   - TTTo2L2Nu2B
   - TTWJets
   - TTZJets
   - TTWWJets
   - TTGJets


The Powheg set
-----------------------------------------------

- qq,qg -> 4l + 1 j (Powheg). NLO. Nickname: qq_ZZ.
   - ZZTo4mu
   - ZZTo4e
   - ZZTo2mu2tau
   - ZZTo2e2tau
   - ZZTo2e2mu
   - ZZTo4tau
- qq -> 4l + 2 j (Phantom samples). LO. Nickname: qq_4l2j.
   - ZZTo4eJJ_SMHContinInterf_H125.6
   - ZZTo2e2muJJ_SMHContinInterf_H125.6
   - ZZTo4muJJ_SMHContinInterf_H125.6
- gg -> 4l (MCFM). LO. Nickname: gg_4l.
   - ggTo4e_SMHContinInterf-MCFM67_H125.6
   - ggTo4mu_SMHContinInterf-MCFM67_H125.6
   - ggTo2e2mu_SMHContinInterf-MCFM67_H125.6
- Higgs samples
   - ttH126
- Dibosons
   - WZ
   - WWJets
   - WGToLNuG
- Tribosons
   - WZZJets (or WZZ_aMCatNLO)
   - WWWJets
   - WWZJets
   - ZZZJets
   - WWGJets
- ttbar
   - TTTo2L2Nu2B
   - TTWJets
   - TTZJets
   - TTWWJets
   - TTGJets


Additional remarks
-----------------------------------------------

- gg_4l are at LO, but Passarino calculated the K factor to go to NNLO. It is an inclusive value, that might depend on jets multiplicity.
The correction could be sizable, as for the inclusive KF it was obtained something of the order of 2 (for a population that account for the 10% of the total).
- An EWK correction for all qq process needs to be applied. Right now, it exist as a KF and the depedence as a function of the number of jets is unknown. Also this correction can be sizable.
- Signal lies in qq_ZZ, qq_4l2j, gg_4l, Higgs samples, WZZJets, ZZZJets.
- Irreducible background in: WWZJets, TTZJets, TTWWJets.
