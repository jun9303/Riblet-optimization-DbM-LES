#!/bin/env python
#-*- coding:utf-8 -*-
#
import os, uuid, sys, socket
import numpy as np

from morphing import *
from lesrun import solverRun

if len(sys.argv[1:]) != 13:
    raise RuntimeError('13 args needed')

testDir = os.getcwd() + '/testspace/' + str(uuid.uuid4().hex) +'/'
os.mkdir(testDir)

print('Current workspace is '+ list(filter(None, testDir.split('/')))[-1])

weights = []
for i in range(1, 11):
    weights.append(float(sys.argv[i]))

# weights = getRandomSamplesOnNSphere(10,1.,1)
wplus = float(sys.argv[11]) # np.random.rand() * 50 # max 50
hplus = float(sys.argv[12]) # np.random.rand() * 50 # max 50
splus = float(sys.argv[13]) # np.random.rand() * 50 # max 50
gridDir = testDir
ibmDir = testDir

morphedShape, weights_norm = rib_dbm(os.getcwd() + '/baseline',weights)

with open(testDir + 'weight.dat', 'w') as f:
    np.savetxt(f, weights_norm, delimiter='\n', fmt='%15.12f')
    f.write(('%15.7f' % wplus) + '\n')
    f.write(('%15.7f' % hplus) + '\n')
    f.write(('%15.7f' % splus))

with open(testDir + 'computed_by_' + socket.gethostname(), 'w') as f:
    f.write('\n')

lesGridPre(morphedShape,wplus,hplus,splus,gridDir)
ibmBodyPre(morphedShape,wplus,hplus,splus,gridDir,ibmDir)
solverRun(gridDir,ibmDir,testDir)
