#!/usr/bin/python

import sys
import re
import matplotlib.pyplot as plt

if len(sys.argv) is not 2:
  print 'Usage:', sys.argv[0], '<parma log file>'
  sys.exit()

parmaLogName = sys.argv[1]
infile = open(parmaLogName)

class base(): pass #HACCCCCCCKKKKK

metrics = []
terms = []
class metric():
  def __init__(self,name):
    self.name = name
    metrics.append(self)
    terms.append(name)

def plot(dat, xlabel, ylabel, name, title):
    for d,s in zip(dat,name):
      print s
      plt.plot(d,label=s)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(title)
    plt.close('all')

def plot2(dat, axis, xlabel, ylabel, name, title):
    color=['b','g','r','c','m','y','b']
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    for d,s,a,c in zip(dat,name,axis,color):
      print s
      if a is 'left':
          ax1.plot(d,color=c,linestyle='--',label=s)
      if a is 'right':
          ax2.plot(d,color=c,label=s)
    ax1.set_xlabel(xlabel)
    ax1.set_xlim(left=-10)
    ax1.set_ylabel(ylabel[0])
    ax1.legend(loc='upper left')
    ax2.set_ylabel(ylabel[1])
    ax2.legend(loc='upper right')
    plt.savefig(title)
    plt.close('all')

vtxImb = metric("vertexImbalance")
vtxImb.has = lambda s: True if s.startswith("vtxImb") or s.startswith("vtx imbalance") else False
vtxImb.get = lambda s: float(s.split()[1]) if s.startswith("vtxImb") else float(s.split()[2])
vtxImb.plot = lambda a,suffix: plot([a],"Iteration","Max Vtx/Avg Vtx", ["vtxImb"], vtxImb.name+suffix)

elmImb = metric("elementImbalance")
elmImb.has = lambda s: True if s.startswith("elmImb") else False
elmImb.get = lambda s: float(s.split()[1])
elmImb.plot = lambda a,suffix: plot(a,"Iteration","Max Elm/Avg Elm", "elmImb", elmImb.name+suffix)

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
bal.plot = lambda a: plot([a],"ParMA Run","Time (seconds)", ["balance"], "balanceTime")

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

rgnsPerPart = runs[0].metrics[rgnAvg.name][0]
print 'average elements per part:', rgnsPerPart

runCount = 0
balAll = []
vtxImbAll = []
elmImbAll = []
iterTimeAll = []
distUpAll = []
elmMigrAll = []
elmSelAll = []
for run in runs:
  for key,value in run.metrics.items():
    if key is vtxImb.name:
      vtxImbAll.extend(value) 
    if key is elmImb.name:
      elmImbAll.extend(value) 
    if key is bal.name:
      balAll.extend(value)
    if key is distUp.name:
      distUpAll.extend(value)
    if key is elmMigr.name:
      elmMigrAll.extend(value)
    if key is elmSel.name:
      elmSelAll.extend(value)
  iterTime = [sum(x) for x in zip(\
                        run.metrics[elmSel.name], \
                        run.metrics[distUp.name], \
                        run.metrics[elmMigr.name])]
  iterTimeAll.extend(iterTime)
  d = len(vtxImbAll) - len(elmImbAll)
  if d > 0:
    elmImbAll.extend([None]*d)
  d = len(vtxImbAll) - len(iterTimeAll)
  if d > 0:
    iterTimeAll.extend([None]*d)
    distUpAll.extend([None]*d)
    elmMigrAll.extend([None]*d)
    elmSelAll.extend([None]*d)
  runCount = runCount + 1

vtxImbPtReduction = (vtxImbAll[0]-vtxImbAll[-1])*100
totTime = sum(balAll)
print 'total number of parma iterations: ' + str(len(vtxImbAll))
print 'percentage point reduction in vertex imbalance: ' + str(vtxImbPtReduction)
print 'total time (s) for parma balancing: ' + str(totTime)
print 'time (s) per parma balance iteration: ' + str(totTime/len(vtxImbAll))
print 'time (s) per vertex imbalance percentage point reduction: ' + str(totTime/vtxImbPtReduction)
print 'vertex imbalance percentage point reduction per element per second: ' + str(vtxImbPtReduction/totTime/rgnsPerPart)

plot2(dat=[iterTimeAll, vtxImbAll, elmImbAll], \
     axis=['right', 'left', 'left'], \
     name=['iterTime', 'vtxImb', 'elmImb'], \
     xlabel="Iteration", \
     ylabel=['Max/Avg', 'Time (s)'], \
     title="vtxElmImbVsIterTime")

bal.plot(balAll)
plot([vtxImbAll, elmImbAll], "Iteration", "Max/Avg", ["vtxImb", "elmImb"], "vtxElmImbalance")
      
infile.close()
