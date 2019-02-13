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
from collections import OrderedDict
from multiprocessing import Pool, Process
#from multiprocessing.dummy import Pool, Process
import pandas, root_numpy
import root_pandas
from itertools import product
gr.SetBatch(True) # NEEDS TO BE SET FOR MULTIPROCESSING OF plot.Draw()
import PhysicsTools.HeppyCore.framework.config as cfg
import os
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
pf.setpfstyle()
####################################################################################################
plotDir     = '/t3home/vstampf/eos/plots/SF/'

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



#####################################################################################################
                                    ##### histos #####                           
#####################################################################################################
l_pt  = [5.0, 10.0, 20.0, 35.0, 50.0, 100.0, 200.0, 500.0]
b_pt  = np.array(l_pt)

l_eta = [0.0, 0.8, 1.444, 1.566, 2.000, 2.500]
b_eta = np.array(l_eta)
l_eta_tag = ['00t08', '08t14', '14t16', '16t20', '20t25']    

b_y    = np.arange(0.,1.,0.1)
framer = rt.TH2F('','',len(b_pt)-1,b_pt,len(b_y)-1,b_y)
framer.GetYaxis().SetRangeUser(0.25, 1.0)
framer.GetXaxis().SetRangeUser(1, 505)
framer.GetXaxis().SetMoreLogLabels()
framer.GetXaxis().SetNoExponent()
framer.SetTitle(';p_{T} [GeV]; Efficiency')
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
# trivial
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
# non-trivial
#####################################################################################################
muonIDs = [# 'isTightMuon' # FIXME not there?!
           'Glb',
           'GlbPT',
           'Medium',
           # 'Medium2016',
           'PF',
           'Tight2012',
           'Loose',
           'Track_HP',
           'VBTF',
           'VBTF_nL8',
           'VBTF_nL9',
           'TM',
           'TMA',
           'TMLST',
           'TMLSAT',
           'TMOST',
           'TMOSL',
           'TMOSTQual',
           'HWWID',]
#####################################################################################################
                                    ##### efficiencies #####                           
#####################################################################################################
def getSF(mode='ele'):

    bch = 16
    batches = int(len(eleIDs)/bch) + 1

    # fill histo's
    for n in range(batches):

        procs = []
        for ID in eleIDs[n*bch:(n+1)*bch]:
            for i,tag in enumerate(l_eta_tag):
                proc = Process(target=fillHistos, args=(mode,ID,tag,i))
                procs.append(proc)
                proc.start()

        for proc in procs:
            proc.join()

#    for ID in eleIDs: 
#        computeEffs(mode,ID,True)
    
    # compute eff's
    for n in range(batches):

        procs = []
        for ID in eleIDs[n*bch:(n+1)*bch]:
            proc = Process(target=computeEffs, args=(mode,ID,True))
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()
#####################################################################################################

#####################################################################################################
def getEff(ID=eleIDs[0], mode='ele'):
    fillHistosOld(mode,ID)
    computeEffs(mode, ID)
#####################################################################################################

#####################################################################################################
def fillHistos(mode,ID,tag,i):

    if mode == 'ele':
        inFileDYtmp = 'DY_MG_EGamma.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tnpEleIDs')
        t = tFile.Get('fitter_tree')

    if mode == 'mu':
        inFileDYtmp = 'DY_MG_Muon.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tpTree')
        t = tFile.Get('fitter_tree')

    eta_cut = '%f < abs(el_sc_eta) & abs(el_sc_eta) < %f'%(l_eta[i],l_eta[i+1])

    cuts_all  =  eta_cut + ' & tag_Ele_pt > 30 & abs(tag_sc_eta) < 2.17 & mcTrue == 1 & abs(mass - 91.19) < 20 & el_q + tag_Ele_q == 0' 
    cuts_all  += ' & el_ecalEnergy * sin( 2 * atan( exp(el_sc_eta) ) )  > 0.5 & abs(el_sc_eta) < 2.5'
    cuts_pass = cuts_all + ' & ' + ID

    # ecalEnergy * sin(superClusterPosition.theta)>5.0 &&  (abs(-log(tan(superClusterPosition.theta/2)))<2.5) 
    # eta = -ln tan theta/2
    # 2 * arctan(exp(- eta)) = theta
    # superClusterPosition.theta = 2 * atan(exp(el_sc_eta)) 

    print '\n\t mode: %s, eta: %s, ID: %s, all entries: %i, passing: %i, avg eff: %.2f' %(mode, tag, ID, t.GetEntries(cuts_all), t.GetEntries(cuts_pass), t.GetEntries(cuts_pass)/t.GetEntries(cuts_all))

    outfile_all = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_all.root' %(mode,ID,tag), 'recreate')
    outfile_all.cd()
    t.Draw( 'el_pt>>ALL(99,5,500)', cuts_all)
    outfile_all.Write()
    outfile_all.Close()

    outfile_pass = rt.TFile.Open(plotDir + 'tmp/pt_eta_%s_%s_%s_pass.root' %(mode,ID,tag), 'recreate')
    outfile_pass.cd()
    t.Draw( 'el_pt>>PASS(99,5,500)', cuts_pass)
    outfile_pass.Write()
    outfile_pass.Close()

    print '\n\t filling done... \n'
#####################################################################################################

#####################################################################################################
def fillHistosOld(mode, ID):

    if mode == 'ele':
        inFileDYtmp = 'DY_MG_EGamma.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tnpEleIDs')
        t = tFile.Get('fitter_tree')

    if mode == 'mu':
        inFileDYtmp = 'DY_MG_Muon.root'
        fin = rt.TFile(inFileDYtmp)
        tFile = fin.Get('tpTree')
        t = tFile.Get('fitter_tree')

    for i in range(len(l_eta)-1):
        eta_cut = '%f < abs(tag_Ele_eta) & abs(tag_Ele_eta) < %f'%(l_eta[i],l_eta[i+1])

        cuts_all  = eta_cut + ''
        cuts_pass = cuts_all + ' & ' + ID + ' == 1'

        print '\n\t mode: %s, eta: %s, ID: %s, all entries: %i, passing: %i, avg eff: %.2f' %(mode, l_eta_tag[i], ID, t.GetEntries(cuts_all), t.GetEntries(cuts_pass), t.GetEntries(cuts_pass)/t.GetEntries(cuts_all))

        t.Draw( 'tag_Ele_pt >> pt_eta_%s_all' %l_eta_tag[i], cuts_all)

        t.Draw( 'tag_Ele_pt >> pt_eta_%s_pass'%l_eta_tag[i], cuts_pass)

    print '\n\t filling done... \n'
#####################################################################################################

#####################################################################################################
def computeEffs(mode, ID, fromFile=False): # TODO make this a function

    if fromFile == True:
        f_in = {}
    
        for tag in l_eta_tag:
            f_in['%s_all' %tag] = rt.TFile(plotDir + 'tmp/pt_eta_%s_%s_%s_all.root' %(mode,ID,tag))
            f_in['%s_pass'%tag] = rt.TFile(plotDir + 'tmp/pt_eta_%s_%s_%s_pass.root'%(mode,ID,tag))

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


    eff_pt_eta_00t08 = rt.TEfficiency(h_pt_eta_00t08_pass_rbn, h_pt_eta_00t08_all_rbn)
    eff_pt_eta_00t08.SetMarkerColor(rt.kGreen+2)
    eff_pt_eta_00t08.SetLineColor(rt.kGreen+2)

    eff_pt_eta_08t14 = rt.TEfficiency(h_pt_eta_08t14_pass_rbn, h_pt_eta_08t14_all_rbn)
    eff_pt_eta_08t14.SetMarkerColor(rt.kBlue+2)
    eff_pt_eta_08t14.SetLineColor(rt.kBlue+2)

    eff_pt_eta_14t16 = rt.TEfficiency(h_pt_eta_14t16_pass_rbn, h_pt_eta_14t16_all_rbn)
    eff_pt_eta_14t16.SetMarkerColor(rt.kRed+2)
    eff_pt_eta_14t16.SetLineColor(rt.kRed+2)

    eff_pt_eta_16t20 = rt.TEfficiency(h_pt_eta_16t20_pass_rbn, h_pt_eta_16t20_all_rbn)
    eff_pt_eta_16t20.SetMarkerColor(rt.kYellow+2)
    eff_pt_eta_16t20.SetLineColor(rt.kYellow+2)

    eff_pt_eta_20t25 = rt.TEfficiency(h_pt_eta_20t25_pass_rbn, h_pt_eta_20t25_all_rbn)
    eff_pt_eta_20t25.SetMarkerColor(rt.kCyan+2)
    eff_pt_eta_20t25.SetLineColor(rt.kCyan+2)

    c_eff = rt.TCanvas('eff','eff'); c_eff.cd()
    framer.Draw()
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
    pf.showlumi(re.sub('passing','',ID))
    pf.showlogopreliminary()
    c_eff.SetLogx()
    c_eff.SetGridx(0)
    c_eff.Modified(); c_eff.Update()
    save(c_eff, 'pt_eta', 'SF', mode, ID) 
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
        framer.Draw()
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









