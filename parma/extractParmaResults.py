#!/usr/bin/python

import sys
import re

if len(sys.argv) is not 2:
  print 'Usage:', sys.argv[0], '<parma log file>'
  sys.exit()

parmaLogName = sys.argv[1]
infile = open(parmaLogName)

outFileName = parmaLogName + ".csv"
outfile = open(outFileName, 'w')

class base(): pass #HACCCCCCCKKKKK

metrics = []
terms = []
class metric():
  def __init__(self,name):
    self.name = name
    metrics.append(self)
    terms.append(name)

vtxImb = metric("vertexImbalance")
vtxImb.has = lambda s: True if s.startswith("vtxImb") or s.startswith("vtx imbalance") else False
vtxImb.get = lambda s: float(s.split()[1]) if s.startswith("vtxImb") else float(s.split()[2])

elmImb = metric("elementImbalance")
elmImb.has = lambda s: True if s.startswith("elmImb") else False
elmImb.get = lambda s: float(s.split()[1])

elmSel = metric("elementSelection")
elmSel.has = lambda s: True if s.startswith("elements selected in") else False
elmSel.get = lambda s: float(s.split()[3])

elmMigr = metric("elementMigration")
elmMigr.has = lambda s: True if s.startswith("elements migrated in") else False
elmMigr.get = lambda s: float(s.split()[3])

distUp = metric("distanceUpdate")
distUp.has = lambda s: True if s.startswith("distance updated in") else False
distUp.get = lambda s: float(s.split()[3])

dcMax = metric("disconnectedMax")
dcMax.has = lambda s: True if "disconnected <max avg>" in s else False
dcMax.get = lambda s: int(s.split()[-2])

dcAvg = metric("disconnectedAvg")
dcAvg.has = lambda s: True if "disconnected <max avg>" in s else False
dcAvg.get = lambda s: float(s.split()[-1])

nbMax = metric("neighborsMax")
nbMax.has = lambda s: True if "neighbors <max avg>" in s else False
nbMax.get = lambda s: int(s.split()[-2])

nbAvg = metric("neighborsAvg")
nbAvg.has = lambda s: True if "neighbors <max avg>" in s else False
nbAvg.get = lambda s: float(s.split()[-1])

getTot = lambda s: int(float(s.split()[-4])) #ehhhhh trim off the decimal
getMax = lambda s: int(float(s.split()[-3])) #ehhhhh trim off the decimal
getMin = lambda s: int(float(s.split()[-2])) #ehhhhh trim off the decimal
getAvg = lambda s: float(s.split()[-1])

vtxTot = metric("vertexTotal")
vtxTot.has = lambda s: True if "weighted vtx <tot max" in s else False
vtxTot.get = getTot

vtxMax = metric("vertexMax")
vtxMax.has = vtxTot.has
vtxMax.get = getMax

vtxMin = metric("vertexMin")
vtxMin.has = vtxTot.has
vtxMin.get = getMin

vtxAvg = metric("vertexAvg")
vtxAvg.has = vtxTot.has
vtxAvg.get = getAvg

rgnTot = metric("rgnTotal")
rgnTot.has = lambda s: True if "weighted rgn <tot max" in s else False
rgnTot.get = getTot

rgnMax = metric("rgnMax")
rgnMax.has = rgnTot.has
rgnMax.get = getMax

rgnMin = metric("rgnMin")
rgnMin.has = rgnTot.has
rgnMin.get = getMin

rgnAvg = metric("rgnAvg")
rgnAvg.has = rgnTot.has
rgnAvg.get = getAvg

ownBdryTot = metric("ownedBoundaryVtxTotal")
ownBdryTot.has = lambda s: True if "owned bdry vtx <tot max" in s else False
ownBdryTot.get = getTot

ownBdryMax = metric("ownedBoundaryVtxMax")
ownBdryMax.has = ownBdryTot.has
ownBdryMax.get = getMax

ownBdryMin = metric("ownedBoundaryVtxMin")
ownBdryMin.has = ownBdryTot.has
ownBdryMin.get = getMin

ownBdryAvg = metric("ownedBoundaryVtxAvg")
ownBdryAvg.has = ownBdryTot.has
ownBdryAvg.get = getAvg

sharBdryTot = metric("sharedBoundaryVtxTotal")
sharBdryTot.has = lambda s: True if "shared bdry vtx <tot max" in s else False
sharBdryTot.get = getTot

sharBdryMax = metric("sharedBoundaryVtxMax")
sharBdryMax.has = sharBdryTot.has
sharBdryMax.get = getMax

sharBdryMin = metric("sharedBoundaryVtxMin")
sharBdryMin.has = sharBdryTot.has
sharBdryMin.get = getMin

sharBdryAvg = metric("sharedBoundaryVtxAvg")
sharBdryAvg.has = sharBdryTot.has
sharBdryAvg.get = getAvg

stovMax = metric("sharedSidesToElementsMax")
stovMax.has = lambda s: True if "sharedSidesToElements <max" in s else False
stovMax.get = lambda s: float(s.split()[-3])

stovMin = metric("sharedSidesToElementsMin")
stovMin.has = stovMax.has
stovMin.get = lambda s: float(s.split()[-2])

stovAvg = metric("sharedSidesToElementsAvg")
stovAvg.has = stovMax.has
stovAvg.get = getAvg

sidesAvg = metric("averageSides")
sidesAvg.has = lambda s: True if " avgSides " in s else False
sidesAvg.get = lambda s: float(s.split()[-1])

bal = metric("balancedIn")
bal.has = lambda s: True if "balanced in" in s else False
bal.get = lambda s: float(s.split()[8])
metrics.remove(bal)

def makeRunLog():
  r = base()
  r.metrics = {term:[] for term in terms}
  return r

def getNewRunLog(runs):
  runs.append(makeRunLog())
  return runs[-1]  # the back of the list
 
runs = []
r = getNewRunLog(runs)
for line in infile:
  for m in metrics:
    if m.has(line):
      r.metrics[m.name].append(m.get(line))
  if bal.has(line):
    r.metrics[bal.name].append(bal.get(line))
    r = getNewRunLog(runs)

for run in runs:
  print '---------'
  for key,value in run.metrics.items():
    print key,len(value),value

infile.close()
outfile.close()
