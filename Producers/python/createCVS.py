#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math, csv

from samples_cfi import samplesVVX as samples


#read xsection file
filein  = open('./python/Xsection8TeV_v2.txt','r')

#prepare the csv file
fileoutcsv = open('python/samples_8TeV.csv','w')
fileoutcsv.write("identifier,process,crossSection = -1,totalEvents = -1,luminosity = -1,execute=True,dataset,user,pattern,splitLevel,tune\n")
csvwriter = csv.writer(fileoutcsv) 

#prepare the py file
fileoutpy = open('python/samples_8TeV.py','w')
fileoutpy.write('samples = [')


for i in range(0,len(samples)-1):
    sample = samples[i][0]
    filein.seek(0)
    foundsampleinfile = 0
    for line in filein:
        line =  line.strip()
        if line.startswith("#"):  continue
        found =  line.find(sample)
        if(found>0):
            foundsampleinfile =+ 1
            spline =  line.split(" ")
            foundsample = False
            comment = False
            xsec = -1
            for column in spline:
                if column == sample:
                    foundsample = True
                elif column == '#':
                    comment = True
                elif foundsample and not comment and not column == '' and not column == '\n' and not column == '1' and not column == 'all' and not column == '\t\t':
                    if xsec < 0:
                        xsec = float(column)    
                    else:
                        xsec = xsec * float(column)
            if foundsample:
                # on LXPLUS need to setup python 2.7
                # setenv PATH ${PATH}:/afs/cern.ch/sw/lcg/external/Python/2.7.2/x86_64-slc5-gcc46-opt/bin
                # setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/cern.ch/sw/lcg/external/Python/2.7.2/x86_64-slc5-gcc46-opt/lib
                lineforpy = samples[i] + (round(xsec,10),)
                fileoutpy.write('{},\n'.format(lineforpy))
                # for csv file reshuffle the order
                lineforcsv = [sample,"",round(xsec,10),"","","",samples[i][2],samples[i][1],samples[i][3],samples[i][4],samples[i][5]]
                csvwriter.writerow(lineforcsv)
                if xsec < 0:
                    print "Warning!",sample,"found in xsection file, but without a valid cross section."

    if foundsampleinfile == 0:
        print "{0:s} not found!".format(sample)
        lineforpy = samples[i] + (-1,)
        fileoutpy.write('{},\n'.format(lineforpy))
        lineforcsv = [sample, "",-1,"","","",samples[i][2],samples[i][1],samples[i][3],samples[i][4],samples[i][5]]
        csvwriter.writerow(lineforcsv)

    if foundsampleinfile >1:
        print "More than one instance for {0:s} has been found!".format(sample)

fileoutpy.write(']')                  


