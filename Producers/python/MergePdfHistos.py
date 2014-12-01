#! /usr/bin/env python
import os

samples=["ZZTo4mu","ZZTo4e","ZZTo2e2mu","ZZTo4muJJ_SMHContinInterf","ZZTo4eJJ_SMHContinInterf","ZZTo2e2muJJ_SMHContinInterf","ggTo2e2mu_SMHContinInterf-MCFM67","ggTo4mu_SMHContinInterf-MCFM67","ggTo4e_SMHContinInterf-MCFM67","ZZZJets","ttH126","WH126","ZH126","VBFH126","WZZJets","powheg15H126"]
for s in samples:
    print s
    bashCommand =  "a='hadd "+s+".root '; for file in `ls -d "+s+"_*/`; do a=$a$file'/PDFStudies.root '; done; `echo $a`"
    print '\n'
    os.system(bashCommand)


