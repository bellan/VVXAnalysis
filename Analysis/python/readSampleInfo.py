#! /usr/bin/env python

import sys, os, commands, math


def readSamplesInfo(infoFilePath = 'data/samples_8TeV_2012.csv', indexBy = 'dataset'):
  """
  Loads the sample information database from the given comma-separated-values
  (csv) file.
  """

  import string

  infoFile                  = open(infoFilePath, "r")
  database                  = {}
  defaults                  = {}
  header                    = None
  for line in infoFile:
    line                    = line.strip()
    if not line:              continue
    if line.startswith("#"):  continue

    data                    = line.split(",")

    # Assign info to the entry indexed by the dataset identifier
    if header:            
      if len(data) != len(header):
        raise ValueError, "Inconsistent number of columns in data '" + line + "', expected header = " + str(header)

      info                  = {}
      for (key, value) in zip(header, data):
        if key.startswith("::"):
          if len(value.strip()) == 0:
            value           = []
          else:
            value           = map(string.strip, value.split(";"))
        info[key]           = value

      index                 = info[indexBy]
      if index in database:
        raise ValueError, "Duplicate entries encountered for %d" % index
      del info[indexBy]
      database[index]       = info

    # Read header information (first line in file)
    else:
      header                = []
      for datum in data:
        if not datum:       break
        if "=" in datum:
          (datum, default)  = map(str.strip, datum.split("="))
          defaults[datum]   = default
        header.append(datum)


  if len(database) == 0:
    raise ValueError, "Invalid information file '" + infoFilePath + "', no entries found."
  return (database, defaults)



def readSampleInfo(sample, infoFilePath = 'samples_8TeV_2012.csv', indexBy = 'dataset'):
  db,defaults = readSamplesInfo()

  if sample in db:
    return db[sample]
  else:
    print "Unknown sample", sample
    #sys.exit(2)
    # FIXME!!!!!
    db['crossSection'] = -1
    return db

def crossSection(sample, infoFilePath = 'samples_8TeV_2012.csv', indexBy = 'dataset'):
  return float(readSampleInfo(sample)['crossSection'])
