from __future__ import division
from ROOT import gROOT as gr
import os, platform
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
#import PhysicsTools.HeppyCore.framework.config as cfg
import os
#from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
pf.setpfstyle()
##################################################################################################################################################################################################
eos          = '/eos/user/v/vstampf/'
if platform.platform() == 'Linux-2.6.32-754.3.5.el6.x86_64-x86_64-with-redhat-6.6-Carbon':
   eos       = '/t3home/vstampf/eos/'
plotDir      = eos+'/plots/SF/'
##################################################################################################################################################################################################
                                    ##### histos #####                           
##################################################################################################################################################################################################
l_pt  = [5.0, 10.0, 20.0, 35.0, 50.0, 100.0, 200.0, 500.0]
b_pt  = np.array(l_pt)

### electrons
l_eta_ele = [0.0, 0.8, 1.444, 1.566, 2.000, 2.500]
l_eta_ele = [-2.500, -2.000, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.000, 2.500]
### muons
l_eta_mu = [0.0, 0.9, 1.2, 2.1, 2.4]

b_eta_mu  = np.array(l_eta_mu)
b_eta_ele = np.array(l_eta_ele)
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
##################################################################################################################################################################################################
                                     #####  IDs  #####                           
##################################################################################################################################################################################################
###### MVATight (3/6/19)
###### https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L39
###### 'passingMVATightNew2 == 1'
##################################################################################################################################################################################################
mvaTight = '(( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.2 + (0.032)*(el_pt-10.))) || ( abs(el_eta) < 0.8 && el_pt >= 25' \
           ' && el_MVA94Xnoiso > 0.68) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.1 + (0.025)*(el_pt-10.)))'\
           ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > 0.475) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 &&' \
           ' el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.1 + (0.028)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25' \
           ' && el_MVA94Xnoiso > 0.32))'
##################################################################################################################################################################################################
# https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_passMVAVLooseTightIP2D.py#L97
CutBase_mvaVLooseTightIP2D = ' tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && el_q*tag_Ele_q < 0 && ((( abs(el_eta) < 0.8 && el_pt >=5 && el_pt < 10'\
                             ' && el_MVA94Xnoiso > (0.488)) || ( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.788 + (0.148/15.)*(el_pt-10.)))'\
                             ' || ( abs(el_eta) < 0.8 && el_pt >= 25 && el_MVA94Xnoiso > -0.640) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=5 && el_pt < 10'\
                             ' && el_MVA94Xnoiso > (-0.045)) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25'\
                             ' && el_MVA94Xnoiso > (-0.850 + (0.075/15.)*(el_pt-10.))) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > -0.775)'\
                             ' || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=5 && el_pt < 10 && el_MVA94Xnoiso > (0.176)) || ( abs(el_eta) >= 1.479'\
                             ' && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.81 + (0.077/15.)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5'\
                             ' && el_pt >= 25 && el_MVA94Xnoiso > -0.733))) && passingTightIP2D == 1'
##################################################################################################################################################################################################
# https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtMVATightNewTightIP2D3DIDEmu.py#L98
CutBase_mvaTightIDEmuTightIP2DTightIP3D = 'tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && el_q*tag_Ele_q < 0 && ((( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25' \
                                          ' && el_MVA94Xnoiso > (0.2 + (0.032)*(el_pt-10.))) || ( abs(el_eta) < 0.8 && el_pt >= 25 && el_MVA94Xnoiso > 0.68)' \
                                          ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.1 + (0.025)*(el_pt-10.)))' \
                                          ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > 0.475) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5' \
                                          ' && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.1 + (0.028)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5' \
                                          ' && el_pt >= 25 && el_MVA94Xnoiso > 0.32)) && passingTightIP2D == 1 && passingTightIP3D == 1 && passingIDEmu == 1)' \
                                          ' && tag_Ele_trigMVA > 0.92 && sqrt( 2*mpfMET*tag_Ele_pt*(1-cos(mpfPhi-tag_Ele_phi))) < 45' # low pT add'_cuts only used for multiIso
##################################################################################################################################################################################################
# https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtMVATightNewTightIP2D3DIDEmuConvIHit0.py#L91
CutBase_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits   = 'tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && el_q*tag_Ele_q < 0 && ((( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25'\
                                                            ' && el_MVA94Xnoiso > (0.2 + (0.032)*(el_pt-10.))) || ( abs(el_eta) < 0.8 && el_pt >= 25 && el_MVA94Xnoiso > 0.68)'\
                                                            ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.1 + (0.025)*(el_pt-10.)))'\
                                                            ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > 0.475) || ( abs(el_eta) >= 1.479'\
                                                            ' && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.1 + (0.028)*(el_pt-10.)))'\
                                                            ' || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25 && el_MVA94Xnoiso > 0.32)) && passingTightIP2D == 1'\
                                                            ' && passingTightIP3D == 1 && passingIDEmu == 1 && (el_mHits==0) && passingConvVeto == 1)'
##################################################################################################################################################################################################
### wrt RECO # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Data2018_v1/etc/config/settings_ele_wrtReco.py#L10-L41
eleIDs_reco = {
    'Run2017_CutBasedVetoNoIso94XV1'                   : 'passingCutBasedVetoNoIso94X == 1'     ,
    'Run2017_CutBasedLooseNoIso94XV1'                  : 'passingCutBasedLooseNoIso94X == 1'    ,
    'Run2017_CutBasedMediumNoIso94XV1'                 : 'passingCutBasedMediumNoIso94X == 1'   ,
    'Run2017_CutBasedTightNoIso94XV1'                  : 'passingCutBasedTightNoIso94X == 1'    ,

    'Run2017_CutBasedVetoNoIso94XV2'                   : 'passingCutBasedVetoNoIso94XV2 == 1'   ,
    'Run2017_CutBasedLooseNoIso94XV2'                  : 'passingCutBasedLooseNoIso94XV2 == 1'  ,
    'Run2017_CutBasedMediumNoIso94XV2'                 : 'passingCutBasedMediumNoIso94XV2 == 1' ,
    'Run2017_CutBasedTightNoIso94XV2'                  : 'passingCutBasedTightNoIso94XV2 == 1'  ,

#    'Run2017_MVAVLooseIP2D'                            : 'passingMVAVLoose == 1 && passingTightIP2D == 1',                                                                    #v0         
#    'Run2017_MVAVLooseFOIP2DIDEmu'                     : 'passingMVAVLooseFO == 1 && passingIDEmu == 1 && passingTightIP2D == 1',                                             #v0 
#    'Run2017_MVATightTightIP2D3D'                      : 'passingMVATight == 1 && passingTightIP2D == 1 && passingTightIP3D == 1'                                             #v0
#    'Run2017_MVATightIP2D3DIDEmu'                      : 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1',                       #v0               
#    'Run2017_MVATightTightIP2D3D'                      :  mvaTight['passingMVATightNew2 == 1'] + ' && passingTightIP2D == 1 && passingTightIP3D == 1',                        #v1 added 3/6   
#    'Run2017_MVATightIP2D3DIDEmu'                      :  mvaTight['passingMVATightNew2 == 1'] + ' && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1',   #v1 added 3/6 
#    'Run2017_LeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04' : 'passingLeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04 == 1',                                                              #v0             
#    'Run2017_LeptonMvaMIDEmuTightIP2DSIP3D8miniIso04'  : 'passingLeptonMvaMIDEmuTightIP2DSIP3D8miniIso04 == 1',                                                               #v0            

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L43
    # added 3/7
    'Run2017_MVAVLooseIP2D'                            : '(( abs(el_eta) < 0.8 && el_pt >=5 && el_pt < 10 && el_MVA94Xnoiso > (0.488)) || ( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.788 + (0.148/15.)*(el_pt-10.))) || ( abs(el_eta) < 0.8 && el_pt >= 25 && el_MVA94Xnoiso > -0.640)' \
                                                         ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=5 && el_pt < 10 && el_MVA94Xnoiso > (-0.045)) || ( abs(el_eta) >= 0.8' \
                                                         ' && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.850 + (0.075/15.)*(el_pt-10.))) || ( abs(el_eta) >= 0.8' \
                                                         ' && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > -0.775) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=5' \
                                                         ' && el_pt < 10 && el_MVA94Xnoiso > (0.176)) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.81 + (0.077/15.)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > -0.733)) && passingTightIP2D == 1', 

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L45
    # added 3/7
    'Run2017_MVAVLooseFOIP2DIDEmu'                     : '(( abs(el_eta) < 0.8 && el_pt >=5 && el_pt < 10 && el_MVA94Xnoiso > (-0.135)) || ( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.930 + (0.043/15.)*(el_pt-10.))) || ( abs(el_eta) < 0.8 && el_pt >= 25 && el_MVA94Xnoiso > -0.887)' \
                                                         ' || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=5 && el_pt < 10 && el_MVA94Xnoiso > (-0.417)) || ( abs(el_eta) >= 0.8' \
                                                         ' && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (-0.93 + (0.04/15.)*(el_pt-10.))) || ( abs(el_eta) >= 0.8' \
                                                         ' && abs(el_eta) < 1.479 && el_pt >= 25 && el_MVA94Xnoiso > -0.890) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=5' \
                                                         ' && el_pt < 10 && el_MVA94Xnoiso > (-0.470)) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.942 + (0.032/15.)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > -0.910)) && passingTightIP2D == 1 && passingIDEmu == 1', 

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L40
    # added 3/7
    'Run2017_MVATightTightIP2D3D'                      : '((( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.2 + (0.032)*(el_pt-10.))) || ( abs(el_eta) < 0.8' \
                                                         ' && el_pt >= 25 && el_MVA94Xnoiso > 0.68) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (0.1 + (0.025)*(el_pt-10.))) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > 0.475) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.1 + (0.028)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > 0.32)) && passingTightIP2D == 1 && passingTightIP3D == 1)',

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L41
    # added 3/7
    'Run2017_MVATightIP2D3DIDEmu'                      : '(( abs(el_eta) < 0.8 && el_pt >=10 && el_pt < 25 && el_MVA94Xnoiso > (0.2 + (0.032)*(el_pt-10.))) || ( abs(el_eta) < 0.8' \
                                                         ' && el_pt >= 25 && el_MVA94Xnoiso > 0.68) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (0.1 + (0.025)*(el_pt-10.))) || ( abs(el_eta) >= 0.8 && abs(el_eta) < 1.479 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > 0.475) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >=10 && el_pt < 25' \
                                                         ' && el_MVA94Xnoiso > (-0.1 + (0.028)*(el_pt-10.))) || ( abs(el_eta) >= 1.479 && abs(el_eta) < 2.5 && el_pt >= 25' \
                                                         ' && el_MVA94Xnoiso > 0.32)) && passingTightIP2D == 1 && passingTightIP3D == 1 && passingIDEmu == 1',}

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L50
    # added 3/7
#    'Run2017_LeptonMvaMIDEmuTightIP2DSIP3D8miniIso04' : '((el_MVATTH>0.85) && passingIDEmu == 1 && passingTightIP2D == 1 && (abs(el_sip3d)<8) && passingMini4 == 1)',    # WILL BE RETRAINED

    # https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtReco.py#L51
    # added 3/7
#    'Run2017_LeptonMvaVTIDEmuTightIP2DSIP3D8miniIso04'  : '((el_MVATTH>0.90) && passingIDEmu == 1 && passingTightIP2D == 1 && (abs(el_sip3d)<8) && passingMini4 == 1)',  # WILL BE RETRAINED


    ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose == 1 && passingTightIP2D == 1'
    ### https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_passMVAVLooseTightIP2D.py
#    'Run2017_MVAVLooseTightIP2DMini'                   : 'passingMVAVLooseMini == 1',   #v0
#    'Run2017_MVAVLooseTightIP2DMini2'                  : 'passingMVAVLooseMini2 == 1',  #v0 
#    'Run2017_MVAVLooseTightIP2DMini4'                  : 'passingMVAVLooseMini4 == 1',  #v0

eleIDs_mvaVLooseTightIP2D = {
    # MiniIsolations: https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_passMVAVLooseTightIP2D.py#L22-L24
    # These are with respect to MVAVLoose + IP2D. You can just use "passingMini*" instead of the above long cut string. The MVA cut is redundant since you put the cut in denominator itself.
    # added 3/7
    'Run2017_MVAVLooseTightIP2DMini'                   : 'passingMini == 1',   
    'Run2017_MVAVLooseTightIP2DMini2'                  : 'passingMini2 == 1',   
    'Run2017_MVAVLooseTightIP2DMini4'                  : 'passingMini4 == 1',  
}

eleIDs_mvaTightIDEmuTightIP2DTightIP3D = {
    ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1'
    ### https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtMVATightNewTightIP2D3DIDEmu.py#L19-L26
#    'Run2017_MultiIsoM'                                : 'passingMultiIsoM == 1',
#    'Run2017_MultiIsoT'                                : 'passingMultiIsoT == 1',
#    'Run2017_MultiIsoEmu'                              : 'passingMultiIsoEmu == 1',
    'Run2017_ConvIHit0'                                : 'passingConvVeto == 1 && el_mHits == 0',
    'Run2017_ConvIHit1'                                : 'passingConvVeto == 1 && el_mHits < 2',
#    'Run2017_MultiIsoNew'                              : '( (el_miniIsoAll/el_pt) < 0.09 ) && ( el_ptRatio > 0.85 || el_ptRel > 9.2 )',
    'Run2017_MultiIsoJECv32'                           : '( (el_miniIsoAll/el_pt) < 0.07 ) && ( el_ptRatio > 0.78 || el_ptRel > 8.0 )',
    'Run2017_MultiIsoEmuJECv32'                        : '( (el_miniIsoAll/el_pt) < 0.07 ) && ( el_ptRatio > 0.78 || el_ptRel > 8.0) && passingISOEmu == 1',
}

eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits = {
    ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1' \
    ###                                                                              ' && passingTightIP3D == 1 && passingConvVeto == 1 && el_mHits == 0'
    ### https://github.com/vhegde91/egm_tnp_analysis/blob/LPC_Moriond18_v3.0/etc/config/settings_ele_wrtMVATightNewTightIP2D3DIDEmuConvIHit0.py
    'Run2017_3Qagree'                                  : 'passingCharge == 1',}

eleIDs = dict(eleIDs_reco.items() + eleIDs_mvaVLooseTightIP2D.items() + eleIDs_mvaTightIDEmuTightIP2DTightIP3D.items() + eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits.items())
##################################################################################################################################################################################################
muonIDs_Loose = {
    'miniIso04_LooseId'                      :   'Probe_passMiniIsoL == 1',
    'miniIso02_LooseId'                      :   'Probe_passMiniIsoM == 1',  
    'miniIso01_LooseId'                      :   'Probe_passMiniIsoT == 1',  
    'miniIso005_LooseId'                     :   'Probe_passMiniIsoVT == 1',  
    'MultiIsoL_LooseId'                      :   'Probe_passMultiIsoL == 1',
    'MultiIsoM_LooseId'                      :   'Probe_passMultiIsoM == 1',
    'MultiIsoM_17data_V32JEC_LooseId'        :   'Probe_passMultiIsoM2017v2 == 1',}



muonIDs_Medium = {
    'miniIso04_MediumId'                     :   'Probe_passMiniIsoL == 1',
    'miniIso02_MediumId'                     :   'Probe_passMiniIsoM == 1',
    'miniIso01_MediumId'                     :   'Probe_passMiniIsoT == 1',
    'miniIso005_MediumId'                    :   'Probe_passMiniIsoVT == 1',
    'MultiIsoL_MediumId'                     :   'Probe_passMultiIsoL == 1',
    'MultiIsoM_MediumId'                     :   'Probe_passMultiIsoM == 1',
    'MultiIsoM_17data_V32JEC_MediumId'       :   'Probe_passMultiIsoM2017v2 == 1',}


muonIDs_MediumPrompt = {              
    'miniIso04_MediumPrompt'                 :   'Probe_passMiniIsoL == 1',
    'miniIso02_MediumPrompt'                 :   'Probe_passMiniIsoM == 1',
    'miniIso01_MediumPrompt'                 :   'Probe_passMiniIsoT == 1',
    'miniIso005_MediumPrompt'                :   'Probe_passMiniIsoVT == 1',
    'MultiIsoL_MediumPrompt'                 :   'Probe_passMultiIsoL == 1',
    'MultiIsoM_MediumPrompt'                 :   'Probe_passMultiIsoM == 1',
    'MultiIsoM_17data_V32JEC_MediumPrompt'   :   'Probe_passMultiIsoM2017v2 == 1',
    'LeptonMVAL_MediumPrompt'                :   'Probe_passMVAL == 1 && Probe_passMiniIsoL == 1',
    'LeptonMVAM_MediumPrompt'                :   'Probe_passMVAM == 1 && Probe_passMiniIsoL == 1',
    'LeptonMVAT_MediumPrompt'                :   'Probe_passMVAT == 1 && Probe_passMiniIsoL == 1',
#    'LeptonMVAVT_MediumPrompt'               :   'Probe_passMVAVT == 1 && Probe_passMiniIsoL == 1',
}
muonIDs = dict(muonIDs_Loose.items() + muonIDs_Medium.items() + muonIDs_MediumPrompt.items())
##################################################################################################################################################################################################
                                    ##### efficiencies #####                           
##################################################################################################################################################################################################
def getSF(mode='ele',yr=17):
    
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
            proc = Process(target=computeSFs2D, args=(mode,ID,yr))
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()
##################################################################################################################################################################################################

##################################################################################################################################################################################################
def makeEffs2D(mode='mu',FAST=True,RDF=True,yr=17):
 
    if mode == 'mu':  b_eta = b_eta_mu 
    if mode == 'ele': b_eta = b_eta_ele 
    cuts_all  = None
    cuts_pass = None
    eta_cut   = None
    fast = 'FullSim_'
    if FAST: 
       fast = 'FastSim_'

    print'\n\tFastSim:', FAST

    if mode == 'ele':
        if FAST == False: MODE = 'full'
        if FAST == True:  MODE = 'fast'
        inDir = eos+'/ntuples/scalefactors/TnP_EGamma_trees_%s_%s/' %(yr,MODE)
        files = glob(inDir+'*.root')
        t = rt.TChain('Events')
        for f in files:
            t.Add(f)
        print '\n\tentries:', t.GetEntries()
        tFile = fin.Get('tnpEleIDs')
        t = tFile.Get('fitter_tree')
        IDs = eleIDs
 
        cuts_all =  'tag_Ele_pt > 30 && abs(tag_sc_eta) < 2.17 && mcTrue == 1 && abs(mass - 91.19) < 30 && el_q * tag_Ele_q < 0' 
        cuts_all += ' && el_ecalEnergy * sin( 2 * atan( exp(el_sc_eta) ) )  > 0.5 && abs(el_sc_eta) < 2.5'
#       this is now used only for the multi-iso ID's
#        if RDF == False:
#            cuts_all += ' && tag_Ele_trigMVA > 0.92 && sqrt( 2*mpfMET*tag_Ele_pt*(1-cos(mpfPhi-tag_Ele_phi))) < 45'

        lep_pt   = 'el_pt'
        #lep_eta  = 'abs_el_sc_eta'
        lep_eta  = 'el_sc_eta'
        
        if RDF == True:
            df = rdf(t)
            #f_all = df.Filter(cuts_all).Define('abs_el_sc_eta', 'abs(el_sc_eta)')

    if mode == 'mu':
        if FAST == False: MODE = 'full'
        if FAST == True:  MODE = 'fast'
        inDir = eos+'/ntuples/scalefactors/TnP_Muon_trees_%s_%s/' %(yr,MODE)
        files = glob(inDir+'*.root')
        t = rt.TChain('Events')
        for f in files:
            t.Add(f)
        print '\n\tentries:', t.GetEntries()
        IDs = muonIDs
 
        if FAST == False:
            cuts_all = 'Probe_isGenMatched == 1 && Probe_charge * Tag_charge < 0' if yr != 18 else 'Probe_charge * Tag_charge < 0'
        if FAST == True:
            cuts_all = 'Probe_charge * Tag_charge < 0'                         
        lep_pt   = 'Probe_pt' 
        lep_eta  = 'abs_probe_eta'

        if RDF == True:
            df = rdf(t)
            f_all = df.Filter(cuts_all).Define('abs_probe_eta', 'abs(Probe_eta)')

    if RDF == True:
        n_bef = df.Count().GetValue()
        n_aft = f_all.Count().GetValue()

    if RDF == False:
        n_bef = t.GetEntries()
        n_aft = t.GetEntries(cuts_all)

    print '\n\tentries before selection: %d'    %n_bef

    print '\n\tentries after pre-selection: %d' %n_aft

    print '\n\tcuts_all: %s\n' %cuts_all


    for ID in IDs.keys():#[:1]:
        filtr_all = '1'

        if ID in eleIDs:
            ## special IDs
            ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose == 1 && passingTightIP2D == 1'
            if ID in eleIDs_mvaVLooseTightIP2D:
                filtr_all = IDs['Run2017_MVAVLooseIP2D'] + ' && ' + CutBase_mvaVLooseTightIP2D
                if RDF == True:
                    df_all = f_all.Filter(filtr_all)

            ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1'
            if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3D:
                # df_all = f_all.Filter('passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1')
                # 3/6/19 gio: mvaTight is new
                filtr_all = mvaTight + ' && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1 && ' + CutBase_mvaTightIDEmuTightIP2DTightIP3D
                if RDF == True:
                    df_all = f_all.Filter(filtr_all) 

            # wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight == 1 && passingIDEmu == 1'\
            #                                                                              ' && passingTightIP2D == 1 && passingTightIP3D == 1 && passingConvVeto == 1 && el_mHits == 0'
            if ID in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits:
                # df_all  = f_all.Filter('passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1 && passingConvVeto == 1 && el_mHits == 0')
                # 3/6/19 gio: mvaTight is new
                filtr_all = mvaTight + ' && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1 && passingConvVeto == 1 && el_mHits == 0' \
                                       ' && ' + CutBase_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits  
                if RDF == True:
                    df_all  = f_all.Filter(filtr_all)
                    
            if RDF == True: 
                if ID not in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits \
                   and ID not in eleIDs_mvaTightIDEmuTightIP2DTightIP3DConvVetoMissHits \
                   and ID not in eleIDs_mvaVLooseTightIP2D:
                    df_all  = f_all

        if ID in muonIDs:
            ## special IDs
            ### wrt MVA VLoose ID + TightIP2D, 'passingMVAVLoose == 1 && passingTightIP2D == 1'
            if ID in muonIDs_Loose:
                filtr_all = 'Probe_passL == 1'
                if RDF == True:
                    df_all  = f_all.Filter(filtr_all)
     
            ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D, 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1 && passingTightIP3D == 1'
            if ID in muonIDs_Medium:
                filtr_all = 'Probe_passM == 1'
                if RDF == True:
                    df_all  = f_all.Filter(filtr_all)
     
            ### wrt MVA Tight ID + ID Emu + TightIP2D + TightIP3D + ConvVeto + MissHits = 0, 'passingMVATight == 1 && passingIDEmu == 1 && passingTightIP2D == 1'\
            ### '&& passingTightIP3D == 1 && passingConvVeto == 1 && el_mHits == 0'
            if ID in muonIDs_MediumPrompt:
                filtr_all = 'Probe_passMP == 1'
                if RDF == True:
                    df_all  = f_all.Filter(filtr_all)

            if RDF == True:
                if ID not in muonIDs_Loose and ID not in muonIDs_Medium and ID not in muonIDs_MediumPrompt:
                    df_all  = f_all

        if RDF == True:
            df_pass = df_all.Filter(IDs[ID])

            if mode == 'mu': 
                h_all_ptr  = df_all .Histo2D(('eta_all', 'eta_all', len(b_pt)-1,b_pt,len(b_eta)-1,b_eta), lep_pt, lep_eta) 
                h_pass_ptr = df_pass.Histo2D(('eta_pass','eta_pass',len(b_pt)-1,b_pt,len(b_eta)-1,b_eta), lep_pt, lep_eta) 

            if mode == 'ele': 
                h_all_ptr  = df_all .Histo2D(('eta_all', 'eta_all' ,len(b_eta)-1,b_eta,len(b_pt)-1,b_pt), lep_eta, lep_pt) 
                h_pass_ptr = df_pass.Histo2D(('eta_pass','eta_pass',len(b_eta)-1,b_eta,len(b_pt)-1,b_pt), lep_eta, lep_pt) 

            h_all  = h_all_ptr.GetPtr()
            h_pass = h_pass_ptr.GetPtr()

        if RDF == False:
            filtr_pass = filtr_all + ' && ' + IDs[ID]

            if mode == 'mu': 
                h_all  = rt.TH2F('eta_all', 'eta_all', len(b_pt)-1,b_pt,len(b_eta)-1,b_eta)
                h_pass = rt.TH2F('eta_pass','eta_pass',len(b_pt)-1,b_pt,len(b_eta)-1,b_eta)

            if mode == 'ele': 
                h_all  = rt.TH2F('eta_all', 'eta_all', len(b_eta)-1,b_eta,len(b_pt)-1,b_pt)
                h_pass = rt.TH2F('eta_pass','eta_pass',len(b_eta)-1,b_eta,len(b_pt)-1,b_pt)

#            t.Draw( 'abs(el_sc_eta):el_pt>>eta_all' , '( ' + cuts_all + ' && ' + filtr_all  + ' ) * ( 2 * (truePU>=20) + 1 * (truePU<20) )' )
#            t.Draw( 'abs(el_sc_eta):el_pt>>eta_pass', '( ' + cuts_all + ' && ' + filtr_pass + ' ) * ( 2 * (truePU>=20) + 1 * (truePU<20) )' )

            if yr == 17:
                CUTS_ALL  = cuts_all + ' && ' + filtr_all 
                CUTS_PASS = cuts_all + ' && ' + filtr_pass

            if yr == 18:
                CUTS_ALL  = (cuts_all + ' && ' + filtr_all ).replace('el_MVA94Xnoiso','el_noIsoMVA94X')
                CUTS_PASS = (cuts_all + ' && ' + filtr_pass).replace('el_MVA94Xnoiso','el_noIsoMVA94X')

#            t.Draw( 'abs(el_sc_eta):el_pt>>eta_all' , CUTS_ALL)
#            t.Draw( 'abs(el_sc_eta):el_pt>>eta_pass', CUTS_PASS)

            t.Draw( 'el_pt:el_sc_eta>>eta_all' , CUTS_ALL)
            t.Draw( 'el_pt:el_sc_eta>>eta_pass', CUTS_PASS)

        if yr == 18:
            ID = ID.replace('2017','2018')

        print '\n\tmode: %s, ID: %s, all entries: %i, passing: %i, avg eff: %.2f' %(mode, ID, h_all.GetEntries(), h_pass.GetEntries(), h_pass.GetEntries()/h_all.GetEntries())

        print '\n\tfiltr_all: %s' %filtr_all

        print '\n\tfiltr_pass: %s\n' %(filtr_all + ' && ' + IDs[re.sub('Run201._','Run2017_',ID)])

        h_pass.Divide(h_all)

        c_eff = rt.TCanvas('eff','eff'); c_eff.cd()
        if mode == 'mu': 
            c_eff.SetLogx()
            h_pass.SetTitle(';Lepton p_{T} [GeV]; |#eta|; Efficiency')
        if mode == 'ele': 
            c_eff.SetLogy()
            h_pass.SetTitle(';super-cluster-#eta; Lepton p_{T} [GeV]; Efficiency')
        h_pass.Draw('colztextE')
        h_pass.SetAxisRange(0.,1.,'Z')
        if mode == 'mu':  
            h_pass.GetXaxis().SetNoExponent()
            h_pass.GetXaxis().SetMoreLogLabels()
        if mode == 'ele':  
            h_pass.GetYaxis().SetNoExponent()
            h_pass.GetYaxis().SetMoreLogLabels()
#            pf.showlogo('CMS')
        pf.showlumi(re.sub('Run201._','',ID + '_FastSim' if FAST else ID + '_FullSim'))
    #    c_eff.SetGridx(0)
        c_eff.Modified(); c_eff.Update()
        save(c_eff, 'pt_eta', 'SUSY', mode, ID + '_FastSim' if FAST else ID + '_FullSim') 
##################################################################################################################################################################################################

##################################################################################################################################################################################################
def computeSFs2D(mode, ID, yr): 
    
    f_in_fast = rt.TFile(plotDir + 'SUSY_%s_eff_pt_eta_%s_FastSim.root'%(mode,ID)).Get('eff')

    if yr == 18:
        ID = ID.replace('2017','2018')

    f_in_full = rt.TFile(plotDir + 'SUSY_%s_eff_pt_eta_%s_FullSim.root'%(mode,ID)).Get('eff')

    h_in_full = f_in_full.GetPrimitive('eta_pass')
    h_in_fast = f_in_fast.GetPrimitive('eta_pass')

    h_in_full.Divide(h_in_fast)

    c_sf = rt.TCanvas('sf','sf'); c_sf.cd()
    if mode == 'mu':
        h_in_full.SetTitle(';p_{T} [GeV]; |#eta|; FullSim/FastSim')
        c_sf.SetLogx()
    if mode == 'ele':
        h_in_full.SetTitle(';super-cluster-#eta; p_{T} [GeV]; FullSim/FastSim')
        c_sf.SetLogy()
    h_in_full.SetAxisRange(0.8,1.2,'Z')
    h_in_full.Draw('colztextE')
    pf.showlumi('SF_'+re.sub('Run201._','',ID))
#    pf.showlogo('CMS')
#    c_sf.SetGridx(0)
    c_sf.Modified(); c_sf.Update()
    save(c_sf, 'pt_eta', 'SUSY', mode, ID) 
##################################################################################################################################################################################################
                                    #####  IDs  #####                           
##################################################################################################################################################################################################
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
##################################################################################################################################################################################################

##################################################################################################################################################################################################
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
##################################################################################################################################################################################################
                                    ##### utilities #####                           
##################################################################################################################################################################################################
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
##################################################################################################################################################################################################

##################################################################################################################################################################################################
