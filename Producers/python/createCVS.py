#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math, csv

from samples_cfi import samplesVVX as samples


#read xsection file
filein  = open('./python/Xsection8TeV_v2.txt','r')

#prepare the csv file
fileoutcsv = open('../TreeAnalysis/data/samples_8TeV.csv','w')
fileoutcsv.write("identifier,crossSection = -99.99,totalEvents = -999,luminosity = -99.99,dataset\n")
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
            newline = []
            for column in spline:
                if column == sample:
                    foundsample = True
                    for j in range(0,6):
                        newline.append(samples[i][j])
                elif column == '#':
                    comment = True
                elif foundsample and not comment and not column == '' and not column == '\n' and not column == '1' and not column == 'all' and not column == '\t\t':
                    newline.append(float(column))
            if len(newline) == 8:
                newline[6] = round(newline[6] * newline.pop(7),10)
            if not len(newline) == 0:
                #print newline
                fileoutpy.write('{},\n'.format(tuple(newline)))
                # for csv file, simplify the output
                lineforcsv = [newline[0],newline[6],"","",newline[2]]
                csvwriter.writerow(lineforcsv)

    if foundsampleinfile == 0:
        print "{0:s} not found!".format(sample)
        fileoutcsv.write("{0:s},,,,{1:s}\n".format(sample, samples[i][2]))
        newline = samples[i] + (-1,)
        fileoutpy.write('{},\n'.format(newline))
    if foundsampleinfile >1:
        print "More than one instance for {0:s} has been found!".format(sample)

fileoutpy.write(']')                  


