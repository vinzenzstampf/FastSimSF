from __future__ import division
from ROOT import gROOT as gr
import os
import re
import ROOT as rt
import numpy as np
import plotfactory as pf
from glob import glob
import sys
from pdb import set_trace
from copy import deepcopy
from os.path import normpath, basename, split
from collections import OrderedDict, Counter
from multiprocessing import Pool, Process
#from multiprocessing.dummy import Pool, Process
from ROOT import RDataFrame as rdf
import pandas, root_numpy
from itertools import product
gr.SetBatch(True) # NEEDS TO BE SET FOR MULTIPROCESSING OF plot.Draw()
import PhysicsTools.HeppyCore.framework.config as cfg
import os
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
pf.setpfstyle()
####################################################################################################
plotDir     = '/t3home/vstampf/eos/plots/SF/'
#####################################################################################################
                                    ##### histos #####                           
#####################################################################################################
l_pt  = [5.0, 10.0, 20.0, 35.0, 50.0, 100.0, 200.0, 500.0]
b_pt  = np.array(l_pt)

l_eta = [0.0, 0.8, 1.444, 1.566, 2.000, 2.500]
b_eta = np.array(l_eta)
l_eta_tag = ['00t08', '08t14', '14t16', '16t20', '20t25']    

b_y1    = np.arange(0.,1.,0.1)
b_y2    = np.arange(0.,2.5,0.1)

eff_framer = rt.TH2F('','',len(b_pt)-1,b_pt,len(b_y1)-1,b_y1)
eff_framer.GetYaxis().SetRangeUser(0.25, 1.0)
eff_framer.GetXaxis().SetRangeUser(1, 505)
eff_framer.GetXaxis().SetMoreLogLabels()
eff_framer.GetXaxis().SetNoExponent()
eff_framer.SetTitle(';p_{T} [GeV]; Efficiency')

sf_framer = rt.TH2F('','',len(b_pt)-1,b_pt,len(b_y2)-1,b_y2)
sf_framer.GetYaxis().SetRangeUser(0.2, 1.8)
sf_framer.GetXaxis().SetRangeUser(1, 505)
sf_framer.GetXaxis().SetMoreLogLabels()
sf_framer.GetXaxis().SetNoExponent()
sf_framer.SetTitle(';p_{T} [GeV]; FastSim / FullSim')
#####################################################################################################

#####################################################################################################
h_pt_eta_00t08_pass = rt.TH1F('pt_eta_00t08_pass','pt_eta_00t08_pass',len(b_pt)-1,b_pt)
h_pt_eta_00t08_all  = rt.TH1F('pt_eta_00t08_all' ,'pt_eta_00t08_all' ,len(b_pt)-1,b_pt)

h_pt_eta_08t14_pass = rt.TH1F('pt_eta_08t14_pass','pt_eta_08t14_pass',len(b_pt)-1,b_pt)
h_pt_eta_08t14_all  = rt.TH1F('pt_eta_08t14_all' ,'pt_eta_08t14_all' ,len(b_pt)-1,b_pt)

h_pt_eta_14t16_pass = rt.TH1F('pt_eta_14t16_pass','pt_eta_14t16_pass',len(b_pt)-1,b_pt)
h_pt_eta_14t16_all  = rt.TH1F('pt_eta_14t16_all' ,'pt_eta_14t16_all' ,len(b_pt)-1,b_pt)

h_pt_eta_16t20_pass = rt.TH1F('pt_eta_16t20_pass','pt_eta_16t20_pass',len(b_pt)-1,b_pt)
h_pt_eta_16t20_all  = rt.TH1F('pt_eta_16t20_all' ,'pt_eta_16t20_all' ,len(b_pt)-1,b_pt)

h_pt_eta_20t25_pass = rt.TH1F('pt_eta_20t25_pass','pt_eta_20t25_pass',len(b_pt)-1,b_pt)
h_pt_eta_20t25_all  = rt.TH1F('pt_eta_20t25_all' ,'pt_eta_20t25_all' ,len(b_pt)-1,b_pt)
#####################################################################################################
                                     #####  IDs  #####                           
#####################################################################################################
### wrt RECO # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Data2018_v1/etc/config/settings_ele_wrtReco.py#L10-L41
eleIDs_reco = {
    'Run2017_CutBasedVetoNoIso94XV1'                   : 'passingCutBasedVetoNoIso94X'     ,
    'Run2017_CutBasedLooseNoIso94XV1'                  : 'passingCutBasedLooseNoIso94X'    ,
    'Run2017_CutBasedMediumNoIso94XV1'                 : 'passingCutBasedMediumNoIso94X'   ,
    'Run2017_CutBasedTightNoIso94XV1'                  : 'passingCutBasedTightNoIso94X'    ,

    'Run2017_CutBasedVetoNoIso94XV2'                   : 'passingCutBasedVetoNoIso94XV2'   ,
    'Run2017_CutBasedLooseNoIso94XV2'                  : 'passingCutBasedLooseNoIso94XV2'  ,
    'Run2017_CutBasedMediumNoIso94XV2'                 : 'passingCutBasedMediumNoIso94XV2' ,
    'Run2017_CutBasedTightNoIso94XV2'                  : 'passingCutBasedTightNoIso94XV2'  ,

    'Run2017_MVAVLooseIP2D'                            : 'passingMVAVLoose & passingTightIP2D',                                      
    'Run2017_MVAVLooseFOIP2DIDEmu'                     : 'passingMVAVLooseFO & passingIDEmu & passingTightIP2D',                           
    'Run2017_MVATightTightIP2D3D'                      : 'passingMVATight & passingTightIP2D & passingTightIP3D',                           
    'Run2017_MVATightIP2D3DIDEmu'                      : 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D',                 
    'Run2017_LeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04' : 'passingLeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04',
    'Run2017_LeptonMvaMIDEmuTightIP2DSIP3D8miniIso04'  : 'passingLeptonMvaMIDEmuTightIP2DSIP3D8miniIso04',}

### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose & passingTightIP2D'
eleIDs_mvaVLooseTightIP2D = {
    'Run2017_MVAVLooseTightIP2DMini'                   : 'passingMVAVLooseMini',
    'Run2017_MVAVLooseTightIP2DMini2'                  : 'passingMVAVLooseMini2',
    'Run2017_MVAVLooseTightIP2DMini4'                  : 'passingMVAVLooseMini4',}

### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'
eleIDs_mvaTightIDEmuTightIP2DTightIP3D = {
    'Run2017_MultiIsoM'                                : 'passingMultiIsoM',
    'Run2017_MultiIsoT'                                : 'passingMultiIsoT',
    'Run2017_MultiIsoEmu'                              : 'passingMultiIsoEmu',
    'Run2017_ConvIHit0'                                : 'passingConvVeto & el_mHits == 0',
    'Run2017_ConvIHit1'                                : 'passingConvVeto & el_mHits < 2',
   #following taken from here, https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtMVATightNewTightIP2D3DIDEmu.py#L23-L26
   # TODO CHECK WITH GIOVANNI
    'Run2017_MultiIsoNew'                              : '( (el_miniIsoAll/el_pt) < 0.09 ) & ( el_ptRatio > 0.85 | el_ptRel > 9.2 )',
    'Run2017_MultiIsoJECv32'                           : '( (el_miniIsoAll/el_pt) < 0.07 ) & ( el_ptRatio > 0.78 | el_ptRel > 8.0 )',
    'Run2017_MultiIsoEmuJECv32'                        : '( (el_miniIsoAll/el_pt) < 0.07 ) & ( el_ptRatio > 0.78 | el_ptRel > 8.0) & passingISOEmu',}

# wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'
eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits = {
    'Run2017_3Qagree'                                  : 'passingCharge',}

eleIDs = dict(eleIDs_reco.items() + eleIDs_mvaVLooseTightIP2D.items() + eleIDs_mvaTightIDEmuTightIP2DTightIP3D.items() + eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits.items())
#####################################################################################################
muonIDs_Loose = {
    'miniIso04_LooseId'         :   'Probe_passMiniIsoL',
    'miniIso02_LooseId'         :   'Probe_passMiniIsoM',  
    'MultiIsoL_LooseId'         :   'Probe_passMultiIsoL',
    'MultiIsoM_LooseId'         :   'Probe_passMultiIsoM',}

muonIDs_Medium = {
    'miniIso04_MediumId'        :   'Probe_passMiniIsoL',
    'miniIso02_MediumId'        :   'Probe_passMiniIsoM',
    'MultiIsoL_MediumId'        :   'Probe_passMultiIsoL',
    'MultiIsoM_MediumId'        :   'Probe_passMultiIsoM',}

muonIDs_MediumPrompt = {
    'miniIso04_MediumPrompt'    :   'Probe_passMiniIsoL',
    'miniIso02_MediumPrompt'    :   'Probe_passMiniIsoL',
    'MultiIsoL_MediumPrompt'    :   'Probe_passMultiIsoL',
    'MultiIsoM_MediumPrompt'    :   'Probe_passMultiIsoM',
    'LeptonMVAL_MediumPrompt'   :   'Probe_passMVAL  & Probe_passMiniIsoL',
    'LeptonMVAM_MediumPrompt'   :   'Probe_passMVAM  & Probe_passMiniIsoL',
    'LeptonMVAT_MediumPrompt'   :   'Probe_passMVAT  & Probe_passMiniIsoL',
#    'LeptonMVAVT_MediumPrompt'  :   'Probe_passMVAVT & Probe_passMiniIsoL',
}
muonIDs = dict(muonIDs_Loose.items() + muonIDs_Medium.items() + muonIDs_MediumPrompt.items())
#####################################################################################################
                                    ##### efficiencies #####                           
#####################################################################################################
def getEff(mode='ele',FAST=False, makeHistos=True):
    
    print 'fast:', FAST

    if mode == 'ele':
        IDs = eleIDs
    if mode == 'mu':
        IDs = muonIDs

    bch = 16
    batches = int(len(IDs)/bch) + 1
 
    # fill histo's
    if makeHistos == True:
        for n in range(batches):

            print '\n\tbatch:', n

            procs = []
            for ID in IDs.keys()[n*bch:(n+1)*bch]:
                for i,tag in enumerate(l_eta_tag):
                    proc = Process(target=fillHistos1D, args=(mode,ID,tag,i,FAST))
                    procs.append(proc)
                    proc.start()

            for proc in procs:
                proc.join()
    
    # compute eff's
    for n in range(batches):

        procs = []
        for ID in IDs.keys()[n*bch:(n+1)*bch]:
            proc = Process(target=computeEffs, args=(mode,ID,FAST,False))
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()
#####################################################################################################

#####################################################################################################
def getSF(mode='ele'):
    
    if mode == 'ele':
        IDs = eleIDs
    if mode == 'mu':
        IDs = muonIDs

    bch = 16
    batches = int(len(IDs)/bch) + 1

    # compute SF's
    for n in range(batches):

        procs = []
        for ID in IDs.keys()[n*bch:(n+1)*bch]:
            proc = Process(target=computeSFs, args=(mode,ID))
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()
#####################################################################################################

#####################################################################################################
def fillHistos1D(mode,ID,tag,i,FAST=False):
 
    cuts_all  = None
    cuts_pass = None
    eta_cut   = None
    fast = ''
    if FAST: fast = 'FS_'

    if mode == 'ele':
        inFileDYtmp = 'DY_MG_EGamma_FS.root' if FAST else 'DY_MG_EGamma.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tnpEleIDs')
        t = tFile.Get('fitter_tree')

        eta_cut       = '%f < abs(el_sc_eta) & abs(el_sc_eta) < %f'%(l_eta[i],l_eta[i+1])
        cuts_all_ele  =  eta_cut + ' & tag_Ele_pt > 30 & abs(tag_sc_eta) < 2.17 & mcTrue == 1 & abs(mass - 91.19) < 20 & el_q * tag_Ele_q < 0' 
        cuts_all_ele  += ' & el_ecalEnergy * sin( 2 * atan( exp(el_sc_eta) ) )  > 0.5 & abs(el_sc_eta) < 2.5'

        # ecalEnergy * sin(superClusterPosition.theta)>5.0 &&  (abs(-log(tan(superClusterPosition.theta/2)))<2.5) 
        # eta = -ln tan theta/2
        # 2 * arctan(exp(- eta)) = theta
        # superClusterPosition.theta = 2 * atan(exp(el_sc_eta)) 

        ## special IDs
        ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose & passingTightIP2D'
        if ID in eleIDs_mvaVLooseTightIP2D:
            cuts_all_ele += ' & passingMVAVLoose & passingTightIP2D'

        ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'
        if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3D:
            cuts_all_ele += ' & passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'

        # wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'
        if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits:
            cuts_all_ele += ' & passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'

        cuts_pass_ele = cuts_all_ele + ' & ' + eleIDs[ID]

        cut_wghd_all  = '( ' + cuts_all_ele  + ' ) * ( 2 * (truePU>=20) + 1 * (truePU<20) )' 
        cut_wghd_pass = '( ' + cuts_pass_ele + ' ) * ( 2 * (truePU>=20) + 1 * (truePU<20) )'
     
        cuts_all  = cuts_all_ele
        cuts_pass = cuts_pass_ele
        lep_pt    = 'el_pt'
    
    if mode == 'mu':
        inFileDYtmp = 'DY_MG_Muon_FS.root' if FAST else 'DY_IDK_MUON_IDK.root'
        fin = rt.TFile(inFileDYtmp)
#        tFile = fin.Get('tpTree')
#        t = tFile.Get('fitter_tree')
        t = fin.Get('Events')

        eta_cut = '%f < abs(Probe_eta) & abs(Probe_eta) < %f'%(l_eta[i],l_eta[i+1])
        cuts_all_mu  = eta_cut + ' & Probe_isGenMatched & Probe_charge * Tag_charge < 0'
        cuts_pass_mu = cuts_all_mu + ' & ' + muonIDs[ID]

        ## special IDs
        ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose & passingTightIP2D'
        if ID in muonIDs_Loose:
            cuts_all_mu += ' & Probe_passL'

        ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'
        if ID in muonIDs_Medium:
            cuts_all_mu += ' & Probe_passM'

        # wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'
        if ID in muonIDs_MediumPrompt:
            cuts_all_mu += ' & Probe_passMP'

        cuts_all  = cuts_all_mu
        cuts_pass = cuts_pass_mu
        lep_pt      = 'Probe_pt' 

#    cuts_all_wghd  = '( ' + cuts_all  + ' ) * %d' %(1./t.GetEntries())
#    cuts_pass_wghd = '( ' + cuts_pass + ' ) * %d' %(1./t.GetEntries())

    print '\n\t mode: %s, eta: %s, ID: %s, all entries: %i, passing: %i, avg eff: %.2f' %(mode, tag, \
           ID, t.GetEntries(cuts_all), t.GetEntries(cuts_pass), t.GetEntries(cuts_pass)/t.GetEntries(cuts_all))

    outfile_all = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_%sall.root' %(mode,ID,tag,fast), 'recreate')
    outfile_all.cd()
    t.Draw( lep_pt+'>>ALL(99,5,500)', cuts_all)
    outfile_all.Write()
    outfile_all.Close()

    outfile_pass = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_%spass.root' %(mode,ID,tag,fast), 'recreate')
    outfile_pass.cd()
    t.Draw( lep_pt+'>>PASS(99,5,500)', cuts_pass)
    outfile_pass.Write()
    outfile_pass.Close()

#    print '\n\t filling done... \n'
#####################################################################################################

#####################################################################################################
def fillHistos(mode,ID,tag,i,FAST=False):
 
    cuts_all  = None
    cuts_pass = None
    eta_cut   = None
    fast = ''
    if FAST: fast = 'FS_'

    if mode == 'ele':
        inFileDYtmp = 'DY_MG_EGamma_FS.root' if FAST else 'DY_MG_EGamma.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tnpEleIDs')
        t = tFile.Get('fitter_tree')
        df = rdf(t)

#        eta_cut       = '%f < abs(el_sc_eta) & abs(el_sc_eta) < %f'%(l_eta[i],l_eta[i+1])
        cuts_all_ele  =  'tag_Ele_pt > 30 & abs(tag_sc_eta) < 2.17 & mcTrue == 1 & abs(mass - 91.19) < 20 & el_q * tag_Ele_q < 0' 
        cuts_all_ele  += ' & el_ecalEnergy * sin( 2 * atan( exp(el_sc_eta) ) )  > 0.5 & abs(el_sc_eta) < 2.5'

        f_all = df.Filter(cuts_all_ele).Define('abs_el_sc_eta', 'abs(el_sc_eta)')

        for ID in IDs:

        ## special IDs
        ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose & passingTightIP2D'
        if ID in eleIDs_mvaVLooseTightIP2D:
            df_all = f_all.Filter('passingMVAVLoose & passingTightIP2D')

        ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'
        if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3D:
            df_all = f_all.Filter('passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D')

        # wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'
        if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits:
            df_all  = f_all.Filter('passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0')
 
        if ID not in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits and not in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits and not in eleIDs_mvaVLooseTightIP2D:
            df_all  = f_all

        df_pass = df_all.Filter(eleIDs[ID])

        h_all  = df_all .Histo2D(('eta_all', 'eta_all', len(b_pt)-1,b_pt,len(b_eta)-1,b_eta), 'el_pt', 'abs_el_sc_eta') 
        h_pass = df_pass.Histo2D(('eta_pass','eta_pass',len(b_pt)-1,b_pt,len(b_eta)-1,b_eta), 'el_pt', 'abs_el_sc_eta') 

        outfile_all = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%sall.root' %(mode,ID,fast), 'recreate')
        outfile_all.cd()
        h_all.Draw()
        outfile_all.Write()
        outfile_all.Close()

        outfile_pass = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%spass.root' %(mode,ID,fast), 'recreate')
        outfile_pass.cd()
        h_pass.Draw()
        outfile_pass.Write()
        outfile_pass.Close()

    
    if mode == 'mu':
        inFileDYtmp = 'DY_MG_Muon_FS.root' if FAST else 'DY_IDK_MUON_IDK.root'
        fin = rt.TFile(inFileDYtmp)
#        tFile = fin.Get('tpTree')
#        t = tFile.Get('fitter_tree')
        t = fin.Get('Events')

        eta_cut = '%f < abs(Probe_eta) & abs(Probe_eta) < %f'%(l_eta[i],l_eta[i+1])
        cuts_all_mu  = eta_cut + ' & Probe_isGenMatched & Probe_charge * Tag_charge < 0'
        cuts_pass_mu = cuts_all_mu + ' & ' + muonIDs[ID]

        ## special IDs
        ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose & passingTightIP2D'
        if ID in muonIDs_Loose:
            cuts_all_mu += ' & Probe_passL'

        ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D'
        if ID in muonIDs_Medium:
            cuts_all_mu += ' & Probe_passM'

        # wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight & passingIDEmu & passingTightIP2D & passingTightIP3D & passingConvVeto & el_mHits == 0'
        if ID in muonIDs_MediumPrompt:
            cuts_all_mu += ' & Probe_passMP'

        cuts_all  = cuts_all_mu
        cuts_pass = cuts_pass_mu
        lep_pt      = 'Probe_pt' 

#    cuts_all_wghd  = '( ' + cuts_all  + ' ) * %d' %(1./t.GetEntries())
#    cuts_pass_wghd = '( ' + cuts_pass + ' ) * %d' %(1./t.GetEntries())

    print '\n\t mode: %s, eta: %s, ID: %s, all entries: %i, passing: %i, avg eff: %.2f' %(mode, tag, \
           ID, t.GetEntries(cuts_all), t.GetEntries(cuts_pass), t.GetEntries(cuts_pass)/t.GetEntries(cuts_all))

    outfile_all = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_%sall.root' %(mode,ID,tag,fast), 'recreate')
    outfile_all.cd()
    t.Draw( lep_pt+'>>ALL(99,5,500)', cuts_all)
    outfile_all.Write()
    outfile_all.Close()

    outfile_pass = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_%spass.root' %(mode,ID,tag,fast), 'recreate')
    outfile_pass.cd()
    t.Draw( lep_pt+'>>PASS(99,5,500)', cuts_pass)
    outfile_pass.Write()
    outfile_pass.Close()

#    print '\n\t filling done... \n'
#####################################################################################################

#####################################################################################################
def computeEffs(mode, ID, FAST=False, teff=False): # TODO make this a function

    fast = ''
    if FAST == True: fast = 'FS_'

    f_in = {}

    for tag in l_eta_tag:
        f_in['%s_all' %tag] = rt.TFile(plotDir + 'tmp/pt_eta_%s_%s_%s_%sall.root' %(mode,ID,tag,fast))
        f_in['%s_pass'%tag] = rt.TFile(plotDir + 'tmp/pt_eta_%s_%s_%s_%spass.root'%(mode,ID,tag,fast))

    h_pt_eta_00t08_all  = f_in['00t08_all' ].Get('ALL')
    h_pt_eta_00t08_all_rbn  = h_pt_eta_00t08_all.Rebin(7,'pt_eta_%s_%s_00t08_all_rbn',b_pt)
    h_pt_eta_00t08_pass = f_in['00t08_pass'].Get('PASS')
    h_pt_eta_00t08_pass_rbn  = h_pt_eta_00t08_pass.Rebin(7,'pt_eta_%s_%s_00t08_pass_rbn',b_pt)

    h_pt_eta_08t14_all  = f_in['08t14_all' ].Get('ALL')
    h_pt_eta_08t14_all_rbn  = h_pt_eta_08t14_all.Rebin(7,'pt_eta_%s_%s_08t14_all_rbn',b_pt)
    h_pt_eta_08t14_pass = f_in['08t14_pass'].Get('PASS')
    h_pt_eta_08t14_pass_rbn  = h_pt_eta_08t14_pass.Rebin(7,'pt_eta_%s_%s_08t14_pass_rbn',b_pt)

    h_pt_eta_14t16_all  = f_in['14t16_all' ].Get('ALL')
    h_pt_eta_14t16_all_rbn  = h_pt_eta_14t16_all.Rebin(7,'pt_eta_%s_%s_14t16_all_rbn',b_pt)
    h_pt_eta_14t16_pass = f_in['14t16_pass'].Get('PASS')
    h_pt_eta_14t16_pass_rbn  = h_pt_eta_14t16_pass.Rebin(7,'pt_eta_%s_%s_14t16_pass_rbn',b_pt)

    h_pt_eta_16t20_all  = f_in['16t20_all' ].Get('ALL')
    h_pt_eta_16t20_all_rbn  = h_pt_eta_16t20_all.Rebin(7,'pt_eta_%s_%s_16t20_all_rbn',b_pt)
    h_pt_eta_16t20_pass = f_in['16t20_pass'].Get('PASS')
    h_pt_eta_16t20_pass_rbn  = h_pt_eta_16t20_pass.Rebin(7,'pt_eta_%s_%s_16t20_pass_rbn',b_pt)

    h_pt_eta_20t25_all  = f_in['20t25_all' ].Get('ALL')
    h_pt_eta_20t25_all_rbn  = h_pt_eta_20t25_all.Rebin(7,'pt_eta_%s_%s_20t25_all_rbn',b_pt)
    h_pt_eta_20t25_pass = f_in['20t25_pass'].Get('PASS')
    h_pt_eta_20t25_pass_rbn  = h_pt_eta_20t25_pass.Rebin(7,'pt_eta_%s_%s_20t25_pass_rbn',b_pt)


    if teff == False:
        eff_pt_eta_00t08 = rt.TH1F('eff_pt_eta_00t08', 'eff_pt_eta_00t08', len(b_pt)-1, b_pt)
        eff_pt_eta_08t14 = rt.TH1F('eff_pt_eta_08t14', 'eff_pt_eta_08t14', len(b_pt)-1, b_pt)
        eff_pt_eta_14t16 = rt.TH1F('eff_pt_eta_14t16', 'eff_pt_eta_14t16', len(b_pt)-1, b_pt)
        eff_pt_eta_16t20 = rt.TH1F('eff_pt_eta_16t20', 'eff_pt_eta_16t20', len(b_pt)-1, b_pt)
        eff_pt_eta_20t25 = rt.TH1F('eff_pt_eta_20t25', 'eff_pt_eta_20t25', len(b_pt)-1, b_pt)

        eff_pt_eta_00t08.Divide(h_pt_eta_00t08_pass_rbn, h_pt_eta_00t08_all_rbn)
        eff_pt_eta_08t14.Divide(h_pt_eta_08t14_pass_rbn, h_pt_eta_08t14_all_rbn)
        eff_pt_eta_14t16.Divide(h_pt_eta_14t16_pass_rbn, h_pt_eta_14t16_all_rbn)
        eff_pt_eta_16t20.Divide(h_pt_eta_16t20_pass_rbn, h_pt_eta_16t20_all_rbn)
        eff_pt_eta_20t25.Divide(h_pt_eta_20t25_pass_rbn, h_pt_eta_20t25_all_rbn)

    if teff == True:
        eff_pt_eta_00t08 = rt.TEfficiency(h_pt_eta_00t08_pass_rbn, h_pt_eta_00t08_all_rbn)
        eff_pt_eta_08t14 = rt.TEfficiency(h_pt_eta_08t14_pass_rbn, h_pt_eta_08t14_all_rbn)
        eff_pt_eta_14t16 = rt.TEfficiency(h_pt_eta_14t16_pass_rbn, h_pt_eta_14t16_all_rbn)
        eff_pt_eta_16t20 = rt.TEfficiency(h_pt_eta_16t20_pass_rbn, h_pt_eta_16t20_all_rbn)
        eff_pt_eta_20t25 = rt.TEfficiency(h_pt_eta_20t25_pass_rbn, h_pt_eta_20t25_all_rbn)

    eff_pt_eta_00t08.SetMarkerColor(rt.kGreen+2)
    eff_pt_eta_00t08.SetLineColor(rt.kGreen+2)
 
    eff_pt_eta_08t14.SetMarkerColor(rt.kBlue+2)
    eff_pt_eta_08t14.SetLineColor(rt.kBlue+2)

    eff_pt_eta_14t16.SetMarkerColor(rt.kRed+2)
    eff_pt_eta_14t16.SetLineColor(rt.kRed+2)
 
    eff_pt_eta_16t20.SetMarkerColor(rt.kYellow+2)
    eff_pt_eta_16t20.SetLineColor(rt.kYellow+2)
 
    eff_pt_eta_20t25.SetMarkerColor(rt.kCyan+2)
    eff_pt_eta_20t25.SetLineColor(rt.kCyan+2)

    c_eff = rt.TCanvas('eff','eff'); c_eff.cd()
    eff_framer.Draw()
    eff_pt_eta_00t08.Draw('same')
    eff_pt_eta_08t14.Draw('same')
    eff_pt_eta_14t16.Draw('same')
    eff_pt_eta_16t20.Draw('same')
    eff_pt_eta_20t25.Draw('same')
    leg = rt.TLegend(0.57, 0.18, 0.80, 0.4)
    leg.AddEntry(eff_pt_eta_00t08, '0.000 < |#eta| < 0.800')
    leg.AddEntry(eff_pt_eta_08t14, '0.800 < |#eta| < 1.444')
    leg.AddEntry(eff_pt_eta_14t16, '1.444 < |#eta| < 1.566')
    leg.AddEntry(eff_pt_eta_16t20, '1.566 < |#eta| < 2.000')
    leg.AddEntry(eff_pt_eta_20t25, '2.000 < |#eta| < 2.500')
    leg.Draw()
    pf.showlumi(re.sub('Run2017_','',ID + '_FS' if FAST else ID + ''))
    pf.showlogopreliminary()
    c_eff.SetLogx()
    c_eff.SetGridx(0)
    c_eff.Modified(); c_eff.Update()
    save(c_eff, 'pt_eta', 'SUSY', mode, ID + '_FS' if FAST else ID + '') 
#####################################################################################################

#####################################################################################################
def computeSFs(mode, ID): 

    f_in = {}; h_in = {}
    
    f_in['full'] = rt.TFile(plotDir + 'SUSY_%s_eff_pt_eta_%s.root'   %(mode,ID)).Get('eff')
    f_in['fast'] = rt.TFile(plotDir + 'SUSY_%s_eff_pt_eta_%s_FS.root'%(mode,ID)).Get('eff')

    for tag in l_eta_tag:
        h_in['%s_full'%tag] = f_in['full'].GetPrimitive('eff_pt_eta_%s'%tag)
        h_in['%s_fast'%tag] = f_in['fast'].GetPrimitive('eff_pt_eta_%s'%tag)

    sf_pt_eta_00t08 = rt.TH1F('sf_pt_eta_00t08', 'sf_pt_eta_00t08', len(b_pt)-1, b_pt)
    sf_pt_eta_08t14 = rt.TH1F('sf_pt_eta_08t14', 'sf_pt_eta_08t14', len(b_pt)-1, b_pt)
    sf_pt_eta_14t16 = rt.TH1F('sf_pt_eta_14t16', 'sf_pt_eta_14t16', len(b_pt)-1, b_pt)
    sf_pt_eta_16t20 = rt.TH1F('sf_pt_eta_16t20', 'sf_pt_eta_16t20', len(b_pt)-1, b_pt)
    sf_pt_eta_20t25 = rt.TH1F('sf_pt_eta_20t25', 'sf_pt_eta_20t25', len(b_pt)-1, b_pt)

    sf_pt_eta_00t08.Divide(h_in['00t08_fast'], h_in['00t08_full'])
    sf_pt_eta_00t08.SetMarkerColor(rt.kGreen+2)
    sf_pt_eta_00t08.SetLineColor(rt.kGreen+2)

    sf_pt_eta_08t14.Divide(h_in['08t14_fast'], h_in['08t14_full'])
    sf_pt_eta_08t14.SetMarkerColor(rt.kBlue+2)
    sf_pt_eta_08t14.SetLineColor(rt.kBlue+2)

    sf_pt_eta_14t16.Divide(h_in['14t16_fast'], h_in['14t16_full'])
    sf_pt_eta_14t16.SetMarkerColor(rt.kRed+2)
    sf_pt_eta_14t16.SetLineColor(rt.kRed+2)

    sf_pt_eta_16t20.Divide(h_in['16t20_fast'], h_in['16t20_full'])
    sf_pt_eta_16t20.SetMarkerColor(rt.kYellow+2)
    sf_pt_eta_16t20.SetLineColor(rt.kYellow+2)

    sf_pt_eta_20t25.Divide(h_in['20t25_fast'], h_in['20t25_full'])
    sf_pt_eta_20t25.SetMarkerColor(rt.kCyan+2)
    sf_pt_eta_20t25.SetLineColor(rt.kCyan+2)

    c_sf = rt.TCanvas('sf','sf'); c_sf.cd()
    sf_framer.Draw()
    sf_pt_eta_00t08.Draw('same')
    sf_pt_eta_08t14.Draw('same')
    sf_pt_eta_14t16.Draw('same')
    sf_pt_eta_16t20.Draw('same')
    sf_pt_eta_20t25.Draw('same')
    leg = rt.TLegend(0.57, 0.18, 0.80, 0.4)
    leg.AddEntry(sf_pt_eta_00t08, '0.000 < |#eta| < 0.800')
    leg.AddEntry(sf_pt_eta_08t14, '0.800 < |#eta| < 1.444')
    leg.AddEntry(sf_pt_eta_14t16, '1.444 < |#eta| < 1.566')
    leg.AddEntry(sf_pt_eta_16t20, '1.566 < |#eta| < 2.000')
    leg.AddEntry(sf_pt_eta_20t25, '2.000 < |#eta| < 2.500')
    leg.Draw()
    pf.showlumi('SF_'+re.sub('Run2017_','',ID))
    pf.showlogo('CMS')
    c_sf.SetLogx()
    c_sf.SetGridx(0)
    c_sf.Modified(); c_sf.Update()
    save(c_sf, 'pt_eta', 'SUSY', mode, ID) 
#####################################################################################################
                                    #####  IDs  #####                           
#####################################################################################################
###  SUSY ELECTRON IDS
'''
CutBased Veto ID V1 (no iso)                                    Reco electrons  SF  Run2017_CutBasedVetoNoIso94XV1
CutBased Loose ID V1 (no iso)                                   Reco electrons  SF  Run2017_CutBasedLooseNoIso94XV1
CutBased Medium ID V1 (no iso)                                  Reco electrons  SF  Run2017_CutBasedMediumNoIso94XV1
CutBased Tight ID V1 (no iso)                                   Reco electrons  SF  Run2017_CutBasedTightNoIso94XV1

CutBased Veto ID V2 (no iso)                                    Reco electrons  SF  Run2017_CutBasedVetoNoIso94XV2
CutBased Loose ID V2 (no iso)                                   Reco electrons  SF  Run2017_CutBasedLooseNoIso94XV2
CutBased Medium ID V2 (no iso)                                  Reco electrons  SF  Run2017_CutBasedMediumNoIso94XV2
CutBased Tight ID V2 (no iso)                                   Reco electrons  SF  Run2017_CutBasedTightNoIso94XV2

MVA VLoose ID + TightIP2D                                       Reco electrons  SF  Run2017_MVAVLooseIP2D
MVA VLooseFO ID + ID Emu + TightIP2D                            Reco electrons  SF  Run2017_MVAVLooseFOIP2DIDEmu
MVA Tight ID + TightIP2D + TightIP3D                            Reco electrons  SF  Run2017_MVATightTightIP2D3D
MVA Tight ID + ID Emu + TightIP2D + TightIP3D                   Reco electrons  SF  Run2017_MVATightIP2D3DIDEmu
LeptonMVA VT + IDEmu + TightIP2D + SIP8 + miniIso<0.4 (JECv6)   Reco electrons  SF  Run2017_LeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04
LeptonMVA M + IDEmu + TightIP2D + SIP8 + miniIso<0.4 (JEC v6)   Reco electrons  SF  Run2017_LeptonMvaMIDEmuTightIP2DSIP3D8miniIso04

MiniIso < 0.1                                                   MVA VLoose ID + TightIP2D   SF  Run2017_MVAVLooseTightIP2DMini
MiniIso < 0.2                                                   MVA VLoose ID + TightIP2D   SF  Run2017_MVAVLooseTightIP2DMini2
MiniIso < 0.4                                                   MVA VLoose ID + TightIP2D   SF  Run2017_MVAVLooseTightIP2DMini4

MultiIsoMedium (JECv6)                                          MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoM
MultiIsoTight (JECv6)                                           MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoT
MultiIsoTight + IsoEmu (JECv6)                                  MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoEmu
MultiIsoNew (JECv6)                                             MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoNew
MultiIso (JECv32)                                               MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoJECv32
MultiIso + IsoEmu (JECv32)                                      MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_MultiIsoEmuJECv32
ConvVeto + MissHits = 0                                         MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_ConvIHit0
ConvVeto + MissHits < 2                                         MVA Tight ID + ID Emu + TightIP2D + TightIP3D   SF  Run2017_ConvIHit1

3ChargeAgreement                                                MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0 SF  Run2017_3Qagree
'''
#####################################################################################################

#####################################################################################################
###  SUSY MUON IDS
'''
miniIso02_LooseId              Mini Iso < 0.2                  Loose ID        
miniIso04_LooseId              Mini Iso < 0.4                  Loose ID        
MultiIsoL_LooseId              Multi Iso Loose                 Loose ID        
MultiIsoM_LooseId              Multi Iso Medium                Loose ID        

MultiIsoL_MediumId             Multi Iso Loose                 Medium ID       
MultiIsoM_MediumId             Multi Iso Medium                Medium ID       
miniIso02_MediumId             Mini Iso < 0.2                  Medium ID       
miniIso04_MediumId             Mini Iso < 0.4                  Medium ID       

MultiIsoL_MediumPrompt         Multi Iso Loose                 MediumPrompt ID     
miniIso02_MediumPrompt         Mini Iso < 0.2                  MediumPrompt ID     
miniIso04_MediumPrompt         Mini Iso < 0.4                  MediumPrompt ID     
MultiIsoM_MediumPrompt         Multi Iso Medium                MediumPrompt ID     
LeptonMVAL_MediumPrompt        Lepton MVA L  + miniIso04       MediumPrompt ID     
LeptonMVAM_MediumPrompt        Lepton MVA M  + miniIso04       MediumPrompt ID     
LeptonMVAT_MediumPrompt        Lepton MVA T  + miniIso04       MediumPrompt ID     
LeptonMVAVT_MediumPrompt       Lepton MVA VT + miniIso04       MediumPrompt ID     
'''
#####################################################################################################
                                    ##### utilities #####                           
#####################################################################################################
def save(knvs, lbl_str, sample='', ch='', rmrk=''):
    knvs.GetFrame().SetLineWidth(0)
    knvs.Modified(); knvs.Update()
    if len(rmrk):
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}_{rmrk}.png' .format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str, rmrk=rmrk))
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}_{rmrk}.pdf' .format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str, rmrk=rmrk))
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}_{rmrk}.root'.format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str, rmrk=rmrk))
    else:
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}.png' .format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str))
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}.pdf' .format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str))
        knvs.SaveAs('{Dir}{smpl}_{ch}_{ttl}_{lbl}.root'.format(Dir=plotDir, smpl=sample, ttl=knvs.GetTitle(), ch=ch, lbl=lbl_str))

def th1(name, bins, xtitle=''):
    h = rt.TH1F('h_%s'%name, name, len(bins)-1, bins)
#    h.name = name
    h.SetTitle('%s; %s; Counts'%(name, xtitle))
    return h

def th2(name, binsX, binsY, xtitle='', ytitle=''):
    h = rt.TH2F('h_%s'%name, name, len(binsX)-1, binsX, len(binsY)-1, binsY)
    h.SetTitle('%s; %s; %s'%(name, xtitle, ytitle))
    return h

def fill(tree, hist, var, cut='', opt=''):
    tree, hist, var, cut, opt
#    tree.Draw('{v} >> h_{h}'.format( v=var, h=hist.GetName() ), cut, opt)
    tree.Draw('{v} >> {h}'.format( v=var, h=hist.GetName() ), cut, opt)
    print '\tvar: {v} \n\tcut: {c}'.format(v=var, c=cut)
    print 'entries: ', hist.GetEntries()
    return hist
    
def draw(hist, mode=1, log=False):
    c = rt.TCanvas(hist.GetName(), hist.GetName())
    if mode == 2:
        hist.Draw('colz') 
        if log == True: c.SetLogz()
    if mode == 1:
        hist.Draw('ep') 
        if log == True: c.SetLogy()
    if mode == 'eff':
        eff_framer.Draw()
        hist.Draw('same')
    pf.showlogoprelimsim('CMS')
    # pf.showTitle('iso_cut = 0.%s'%iso_str)
    pf.showTitle(hist.GetName())
    c.Modified; c.Update()
    return c

def plot(tupel, name, var, binsX, binsY=[], xtitle='', ytitle='', mode=1, cut='', log=False, opt='', iso=0.15, eta_bin=['full', '']):
    sample_dir, cutuple = tupel
    eta = eta_bin[0]
    eta_cut = eta_bin[1]
    cut_name = cutuple[0]
    fin = rt.TFile(inDir + sample_dir + suffix)
    t = fin.Get('tree')

    if len(cutuple[1]) > 3: cut += ' & ' + cutuple[1]
    if len(eta_cut) > 3:    cut += ' & ' + eta_cut
    ch     = basename(split(normpath(sample_dir))[0]) 
    sample = basename(normpath(sample_dir))
    if mode == 1: 
        hist = th1(name, binsX, xtitle)
    if mode == 2: 
        hist = th2(name, binsX, binsY, xtitle, ytitle)
    if mode == 'eff':
        numer = th1('%s_n'%name, binsX, xtitle)
        denom = th1('%s_d'%name, binsX, xtitle)
        # TODO finish this mode

    print '\nsample name: {s}_{ch}, entries: {n}'.format(s=sample, ch=ch, n=t.GetEntries())
    print '\tfilling hist: {h}'.format(h=hist.GetName())
    filled_hist = fill(t, hist, var, cut, opt)
    print '\thist: {h} entries: {n}\n'.format(h=hist.GetName(), n=filled_hist.GetEntries())
    c = draw(filled_hist, mode, log)
    save(c, iso, sample, ch, eta)
    return filled_hist, c









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


