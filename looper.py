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

tFile = fin.Get('tnpEleIDs')

t = tFile.Get('fitter_tree')

#####################################################################################################
                                    ##### histos #####                           
#####################################################################################################
l_pt  = [5.0, 10.0, 20.0, 35.0, 50.0, 100.0, 200.0, 500.0]
b_pt  = np.array(l_pt)

l_eta = [0.0, 0.8, 1.444, 1.566, 2.000, 2.500]
b_eta = np.array(l_eta)

h_eta_00t08_pass = rt.TH1F('pt_eta_00t08_pass','pt_eta_00t08_pass',len(b_pt)-1,b_pt)
h_eta_00t08_all  = rt.TH1F('pt_eta_00t08_all ','pt_eta_00t08_all ',len(b_pt)-1,b_pt)

h_eta_08t14_pass = rt.TH1F('pt_eta_08t14_pass','pt_eta_08t14_pass',len(b_pt)-1,b_pt)
h_eta_08t14_all  = rt.TH1F('pt_eta_08t14_all ','pt_eta_08t14_all ',len(b_pt)-1,b_pt)

h_eta_14t16_pass = rt.TH1F('pt_eta_14t16_pass','pt_eta_14t16_pass',len(b_pt)-1,b_pt)
h_eta_14t16_all  = rt.TH1F('pt_eta_14t16_all ','pt_eta_14t16_all ',len(b_pt)-1,b_pt)

h_eta_16t20_pass = rt.TH1F('pt_eta_16t20_pass','pt_eta_16t20_pass',len(b_pt)-1,b_pt)
h_eta_16t20_all  = rt.TH1F('pt_eta_16t20_all ','pt_eta_16t20_all ',len(b_pt)-1,b_pt)

h_eta_20t25_pass = rt.TH1F('pt_eta_20t25_pass','pt_eta_20t25_pass',len(b_pt)-1,b_pt)
h_eta_20t25_all  = rt.TH1F('pt_eta_20t25_all ','pt_eta_20t25_all ',len(b_pt)-1,b_pt)
#####################################################################################################
                                     #####  IDs  #####                           
#####################################################################################################
eleIDs=['passingCharge',
        'passingConvIHit0',
        'passingConvIHit0Chg',
        'passingConvIHit1',
        'passingConvVeto',
        'passingCutBasedLooseNoIso',
        'passingCutBasedLooseNoIso94X',
        'passingCutBasedLooseNoIso94XV2',
        'passingCutBasedMediumMini',
        'passingCutBasedMediumMini94X',
        'passingCutBasedMediumNoIso',
        'passingCutBasedMediumNoIso94X',
        'passingCutBasedMediumNoIso94XV2',
        'passingCutBasedStopsDilepton',
        'passingCutBasedTightMini',
        'passingCutBasedTightMini94X',
        'passingCutBasedTightNoIso',
        'passingCutBasedTightNoIso94X',
        'passingCutBasedTightNoIso94XV2',
        'passingCutBasedVetoNoIso',
        'passingCutBasedVetoNoIso94X',
        'passingCutBasedVetoNoIso94XV2',
        'passingFOID2D',
        'passingHLTsafe',
        'passingIDEmu',
        'passingIHit0',
        'passingIHit1',
        'passingISOEmu',
        'passingLeptonMvaM',
        'passingLeptonMvaMIDEmuTightIP2DSIP3D8miniIso04',
        'passingLeptonMvaVT',
        'passingLeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04',
        'passingLoose2D',
        'passingLoose80X',
        'passingLoose94X',
        'passingLoose94XV2',
        'passingMVA80Xwp80',
        'passingMVA80Xwp90',
        'passingMVA94Xwp80iso',
        'passingMVA94Xwp80isoV2',
        'passingMVA94Xwp80noiso',
        'passingMVA94Xwp80noisoV2',
        'passingMVA94Xwp90iso',
        'passingMVA94Xwp90isoV2',
        'passingMVA94Xwp90noiso',
        'passingMVA94Xwp90noisoV2',
        'passingMVA94XwpHZZisoV2',
        'passingMVA94XwpLiso',
        'passingMVA94XwpLisoV2',
        'passingMVA94XwpLnoiso',
        'passingMVA94XwpLnoisoV2',
        'passingMVATight',
        'passingMVAVLoose',
        'passingMVAVLooseFO',
        'passingMVAVLooseMini',
        'passingMVAVLooseMini2',
        'passingMVAVLooseMini4',
        'passingMVAWP80',
        'passingMVAWP90',
        'passingMedium80X',
        'passingMedium94X',
        'passingMedium94XV2',
        'passingMini',
        'passingMini2',
        'passingMini4',
        'passingMultiIsoEmu',
        'passingMultiIsoM',
        'passingMultiIsoT',
        'passingMultiIsoVT',
        'passingTight2D3D',
        'passingTight80X',
        'passingTight94X',
        'passingTight94XV2',
        'passingTightConvIHit0',
        'passingTightID2D3D',
        'passingTightIP2D',
        'passingTightIP3D',
        'passingVeto80X',
        'passingVeto94X',
        'passingVeto94XV2',]
#####################################################################################################
                                    ##### filling #####                           
#####################################################################################################
cuts      = ''

ID = eleIDs[0]
cuts_pass = cuts + ' & ' + ID

t.Draw( P T '>> pt_eta_00t08_all' , cuts_all)
t.Draw( P T '>> pt_eta_00t08_pass', cuts_pass)
#####################################################################################################
                                    ##### efficiencies #####                           
#####################################################################################################
eff_eta_00t08    = rt.TEfficiency(h_pt_eta_00t08_pass, pt_eta_00t08_all)
eff_eta_08t14    = rt.TEfficiency(h_pt_eta_08t14_pass, pt_eta_08t14_all)
eff_eta_14t16    = rt.TEfficiency(h_pt_eta_14t16_pass, pt_eta_14t16_all)
eff_eta_16t20    = rt.TEfficiency(h_pt_eta_16t20_pass, pt_eta_16t20_all)
eff_eta_20t25    = rt.TEfficiency(h_pt_eta_20t25_pass, pt_eta_20t25_all)





















