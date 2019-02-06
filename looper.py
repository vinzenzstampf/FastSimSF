from __future__ import division
from ROOT import gROOT as gr
import os
import ROOT as rt
import numpy as np
import plotfactory as pf
from glob import glob
import sys
from pdb import set_trace
from copy import deepcopy
from os.path import normpath, basename, split
from collections import OrderedDict
from multiprocessing import Pool
from multiprocessing.dummy import Pool
import pandas, root_numpy
import root_pandas
from itertools import product
#gr.SetBatch(True) # NEEDS TO BE SET FOR MULTIPROCESSING OF plot.Draw()
import PhysicsTools.HeppyCore.framework.config as cfg
import os
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
pf.setpfstyle()
####################################################################################################
plotDir     = '/eos/user/v/vstampf/ntuples/DDE_v2/'

### from Vinay Hegde (EGamma)
inFileDYMG_EG  = 'root://cms-xrd-global.cern.ch//store/user/vhegde/EGamma_ntuples/Run2016_17Jul2018_MiniAODv3_TreeV1/mc/TnPTree_mc_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_allExt.root'
inFileDYAMC_EG = 'root://cms-xrd-global.cern.ch//store/user/vhegde/EGamma_ntuples/Run2016_17Jul2018_MiniAODv3_TreeV1/mc/TnPTree_mc_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_all.root'

creator = ComponentCreator()

### from Lesya 
DY_17 = creator.makeMCComponent(
    name    = 'DY_FS_17', 
    dataset = '/DYJetsToLL_M-50_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_pilot_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM',
    user    = 'CMS', 
    pattern = '.*root', 
    useAAA  = True
)

TT_17 = creator.makeMCComponent(
    name    = 'TT_FS_17', 
    dataset = '/TTJets_DiLept_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_lhe_94X_mc2017_realistic_v15-v1/MINIAODSIM',
    user    = 'CMS', 
    pattern = '.*root', 
    useAAA  = True
)

## 16 = ReMiniAOD
DY_16 = creator.makeMCComponent(
    name    = 'DY_FS_16', 
    dataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUSummer16v3Fast_lhe_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM',
    user    = 'CMS', 
    pattern = '.*root', 
    useAAA  = True
)

TT_16 = creator.makeMCComponent(
    name    = 'TT_FS_16', 
    dataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUSummer16v3Fast_lhe_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM',
    user    = 'CMS', 
    pattern = '.*root', 
    useAAA  = True
)


inFileDYtmpEG = 'DY_MG_FS_EG_NTUPLE.root'

fin = rt.TFile(inFileDYtmpEG)


