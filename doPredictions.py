import ROOT as rt
from array import array
import numpy as np

doPunchThroughPred = True

def getPred(nA, nB, nC):
    relUncA = np.sqrt(yieldA)/yieldA if yieldA>0 else 1.0
    relUncB = np.sqrt(yieldB)/yieldB if yieldB>0 else 1.0
    relUncC = np.sqrt(yieldC)/yieldC if yieldC>0 else 1.0
    
    if nA==0:
        return np.inf, np.inf
    else:
        predD = nB*nC/nA
        nD = predD if predD>0 else 1.0
        uncD = nD*np.sqrt(relUncA*relUncA + relUncB*relUncB + relUncC*relUncC)
        return predD, uncD

def getPredPT(yieldD, predD, predDUnc, y0, y0unc):
    diff = yieldD - predD
    nD = yieldD if yieldD>0 else 1.0
    diffUnc = np.sqrt(nD + predDUnc*predDUnc)

    ratio = y0/(1-y0)
    ratioUnc = 1/pow(1-y0,2)*y0unc
    
    predPT = ratio*diff
    predPTUnc = predPT*np.sqrt(pow(ratioUnc/ratio,2) + pow(diffUnc/diff,2))

    return predPT, predPTUnc

filename = "outData_ABCD_2023SR"
dataFile = rt.TFile.Open(filename+".root","READ")
if(doPunchThroughPred):
    dataFilePT = rt.TFile.Open(filename+"_noJetMET.root","READ")
dataHistos = {}

xbins = array('d', [0,10,40,70,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500])
f1 = rt.TF1("f1","expo(0)+[2]",10,500)
f2 = rt.TF1("f2","[0]",10,500)

f1.SetParLimits(0,-10,-1)
f1.SetParLimits(1,-1,0)
f1.SetParLimits(2,0,0.2)

years = ["2016","2017","2018"]

pref1 = "h_dtRechitClusterSize_dPhiClusterMET_fullSelectionWithNHFLeadCut"
pref2 = "h_dtRechitClusterSize_dPhiClusterMET_fullSelection"

sels = {'invertedMB1_MB2': 'MB2withMB1CR',
        'invertedMB1_MB3': 'MB1HitsCR_MB3',
        'invertedMB1_MB4': 'MB1HitsCR_MB4',
        'SR_MB2': 'MB2CR',
        'SR_MB3': 'SRMB3',
        'SR_MB4': 'SRMB4'}

c1 = rt.TCanvas()

print(filename)
for name,sel in sels.items():
    if sel[:2]=='SR':
        dataHistos[name] = dataFile.Get(pref2+'_'+sel+'_2018')
        dataHistos[name].Add(dataFile.Get(pref2+'_'+sel+'_2017'))
        dataHistos[name].Add(dataFile.Get(pref2+'_'+sel+'_2016'))
    else:
        dataHistos[name] = dataFile.Get(pref1+'_'+sel+'_2018')
        dataHistos[name].Add(dataFile.Get(pref1+'_'+sel+'_2017'))
        dataHistos[name].Add(dataFile.Get(pref1+'_'+sel+'_2016'))

    yieldA = dataHistos[name].Integral(0,10,21,-1)
    yieldB = dataHistos[name].Integral(11,-1,21,-1)
    yieldC = dataHistos[name].Integral(0,10,0,20)
    yieldD = dataHistos[name].Integral(11,-1,0,20)

    predD, uncD = getPred(yieldA, yieldB, yieldC)

    print(name, yieldA, yieldB, yieldC, round(predD,1), round(uncD,1), yieldD)

    if name[:3]=='inv':
        station = name[-3:]
        dataHistos[name+"PT"] = dataFilePT.Get("h_matchedJetPt_nMB1Match_"+station+"_invertedJetVetoLooseMuonNoClusterMET_2018")
        dataHistos[name+"PT"].Add(dataFilePT.Get("h_matchedJetPt_nMB1Match_"+station+"_invertedJetVetoLooseMuonNoClusterMET_2017"))
        dataHistos[name+"PT"].Add(dataFilePT.Get("h_matchedJetPt_nMB1Match_"+station+"_invertedJetVetoLooseMuonNoClusterMET_2016"))

        h_p = dataHistos[name+"PT"].ProjectionX("h_p",2,-1)
        h_t = dataHistos[name+"PT"].ProjectionX("h_t",1,-1)

        h_pp = h_p.Rebin(20,"h_pp",xbins)
        h_tt = h_t.Rebin(20,"h_tt",xbins)

        eff = rt.TEfficiency(h_pp,h_tt)
        if(station=='MB2'):
            ytitle = 'MB1 veto efficiency'
        else:
            ytitle = 'MB1 & MB2 veto efficiency'
        xtitle = 'matched jet pT [GeV]'
        eff.SetTitle(";"+xtitle+";"+ytitle)

        f1.SetParameters(-1.4,-0.01,0.04)
        eff.Fit(f1,"s")
        #if(r.Status()!=0):
        #    r = eff.Fit(f2,"s")
        '''
        if(r.NPar()==3):
            y0 = np.exp(r.Parameter(0)) + r.Parameter(2)
        y0unc = np.sqrt(pow(np.exp(r.Parameter(0))*r.ParError(0),2) + pow(r.ParError(2),2))
        else:
            y0 = r.Parameter(0)
            y0unc = r.ParError(0)
        '''
        y0 = np.exp(f1.GetParameter(0)) + f1.GetParameter(2)
        y0unc = np.sqrt(pow(np.exp(f1.GetParameter(0))*f1.GetParError(0),2) + pow(f1.GetParError(2),2))
        predPT, predPTUnc = getPredPT(yieldD, predD, uncD, y0, y0unc)
        
        print(round(y0,3),round(y0unc,3),round(predPT,1),round(predPTUnc,1))


        '''grint = rt.TGraphErrors(20)
        #grint->SetTitle("Fitted line with .95 conf. band");
        for i in range(21):
            #print(xbins[i])
            grint.SetPoint(i, xbins[i], 0)
        #print(grint)
        #r.GetConfidenceIntervals(f1)
        #hint = rt.TH1D("hint","Fitted Gaussian with .95 conf.band", 20, 0, 500)
        ci = np.zeros_like(xbins)
        #rt.TVirtualFitter.GetFitter().GetConfidenceIntervals(grint,0.68)
        rt.TVirtualFitter.GetFitter().GetConfidenceIntervals(20,1,xbins,ci,0.68)
        '''
        eff.SetFillColorAlpha(rt.kRed,0.2)
        eff.Draw()
        f1.SetLineColor(rt.kBlue)
        f1.Draw("E3 same")
        #r.GetConfidenceIntervals(r.FittedBinData)
        c1.Update()
        c1.WaitPrimitive()
        #c1.SaveAs("vetoEff_"+name+".pdf")
