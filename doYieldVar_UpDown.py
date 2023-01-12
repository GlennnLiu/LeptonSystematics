import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)


lumi = 35.92 # Set lumi to be used for MC scaling

folder = '/eos/user/g/geliu/LepUni/BigTrees/Simp/' # Define input folder name
file_name = '/HTauTauHMuMu_unc.root' # Define input dile name

correlated_leptons = True
uncorrelated_leptons = True
corr_factor = 0

# List of samples to run on
List = [
'ggH125',
#'ZZTo4l',
#'ggTo4l'
]

def sigma_event(rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4):
    rez = (sigma1/SF1)**2 + (sigma2/SF2)**2 + (sigma3/SF3)**2 + (sigma4/SF4)**2 + 2*rho*( sigma1*sigma2/SF1/SF2 + sigma1*sigma3/SF1/SF3 + sigma1*sigma4/SF1/SF4 + sigma2*sigma3/SF2/SF3 + sigma2*sigma4/SF2/SF4 + sigma3*sigma4/SF3/SF4)
    if rez < 0.000001:
        return 0
    else:
        return rez

for Type in List:
    # Open the root file and get tree
    mc=TFile.Open(folder+Type+file_name)
    tree = mc.Get("ZZTree/candTree")
    counters = mc.Get("ZZTree/Counters")
    NGen = counters.GetBinContent(40)

    chans=['2e2mu','2eetau','2emutau','2etautau','2mu2mu','2muetau','2mumutau','2mutautau']
    Z1Flav=[-121,-121,-121,-121,-169,-169,-169,-169]
    Z2Flav=[-169,-165,-195,-225,-169,-165,-195,-225]
    
    nom_yield = []
    br_data = 0
    br_fs = []
    yield_e_up = []
    yield_e_dn = []
    uncor_e_up = []
    uncor_e_dn = []
    yield_mu_up = []
    yield_mu_dn = []
    uncor_mu_up = []
    uncor_mu_dn = []
    for ichan in range(len(chans)):
        nom_yield.append(0.)
        br_fs.append(0.)
        yield_e_up.append([])
        yield_e_dn.append([])
        uncor_e_up.append([])
        uncor_e_dn.append([])
        yield_mu_up.append([])
        yield_mu_dn.append([])
        uncor_mu_up.append([])
        uncor_mu_dn.append([])

    # Uncrrelate effect of trigger/reco/selection
    variation_name = []
    variation_name.append("TRIGGER")
    variation_name.append("RECO+SEL")
    for ichan in range(len(chans)):
        for k in range(2):
            yield_e_up[ichan].append(0.)
            yield_e_dn[ichan].append(0.)
            uncor_e_up[ichan].append(0.)
            uncor_e_dn[ichan].append(0.)
            yield_mu_up[ichan].append(0.)
            yield_mu_dn[ichan].append(0.)
            uncor_mu_up[ichan].append(0.)
            uncor_mu_dn[ichan].append(0.)


    for event in tree:# Loop over all events in tree
        
        br_data+=1

        if (abs(tree.LepLepId[2])==15 and (tree.LepPt[2]<20 or abs(tree.LepEta[2])>2.3)) or (abs(tree.LepLepId[3])==15 and (tree.LepPt[3]<20 or abs(tree.LepEta[3])>2.3)):
            continue
        
        findFlav=False
        for ichan,chan in enumerate(chans):
            if tree.Z1Flav==Z1Flav[ichan] and tree.Z2Flav==Z2Flav[ichan]:
                br_fs[ichan]+=1
                findFlav=True
                break
        if not findFlav:
            continue

        if (br_data % int((tree.GetEntries()/10)) == 0):
            print "{} %".format(str(100*br_data/tree.GetEntries() + 1))


        # Calculate nominal weigh using central value of SF
        weight_nom = event.overallEventWeight*1000*lumi*event.xsec/NGen
        # Nominal value of total SF, product of 4 lepton nominal SF
        SF_tot_nom = event.dataMCWeight
        
        nom_yield[ichan]+=weight_nom

        SF_lep_trig = []
        err_lep_trig_up = []
        err_lep_trig_dn = []
        SF_lep_sel = []
        err_lep_sel_up = []
        err_lep_sel_dn = []

        for i in range (4):
            SF_lep_trig.append(1.)
            err_lep_trig_up.append(0.)
            err_lep_trig_dn.append(0.)
            SF_lep_sel.append(event.LepSF[i])
            err_lep_sel_up.append(event.LepSF_Unc[i])
            err_lep_sel_dn.append(event.LepSF_Unc[i])


        for k in range (2):
        # Calculate each variation independently
            if(k==0):
                TRIG = 1
                SEL  = 0
            elif(k==1):
                TRIG = 0
                SEL  = 1

            SF_var_e_up = 1.
            SF_var_e_dn = 1.
            SF_var_mu_up = 1.
            SF_var_mu_dn = 1.

            # Vary SF of each lepton up and down
            for i in range (4):
                if ( correlated_leptons ):
                    if abs(tree.LepLepId[i])==11:
                        SF_var_e_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep_sel[i] + SEL*err_lep_sel_up[i])
                        SF_var_e_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep_sel[i] - SEL*err_lep_sel_dn[i])
                        SF_var_mu_up *= SF_lep_trig[i] * SF_lep_sel[i]
                        SF_var_mu_dn *= SF_lep_trig[i] * SF_lep_sel[i]
                    elif abs(tree.LepLepId[i])==13:
                        SF_var_mu_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep_sel[i] + SEL*err_lep_sel_up[i])
                        SF_var_mu_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep_sel[i] - SEL*err_lep_sel_dn[i])
                        SF_var_e_up *= SF_lep_trig[i] * SF_lep_sel[i]
                        SF_var_e_dn *= SF_lep_trig[i] * SF_lep_sel[i]
                    else:
                        SF_var_e_up *= SF_lep_trig[i] * SF_lep_sel[i]
                        SF_var_e_dn *= SF_lep_trig[i] * SF_lep_sel[i]
                        SF_var_mu_up *= SF_lep_trig[i] * SF_lep_sel[i]
                        SF_var_mu_dn *= SF_lep_trig[i] * SF_lep_sel[i]
                        
            if (uncorrelated_leptons):
                uncor_e_up[ichan][k] += TRIG*(sigma_event(corr_factor, SF_lep_trig[0], SF_lep_trig[1], SF_lep_trig[2], SF_lep_trig[3], err_lep_trig_up[0] if abs(tree.LepLepId[0])==11 else 0., err_lep_trig_up[1] if abs(tree.LepLepId[1])==11 else 0., err_lep_trig_up[2] if abs(tree.LepLepId[2])==11 else 0., err_lep_trig_up[3] if abs(tree.LepLepId[3])==11 else 0.)) + SEL*(sigma_event(corr_factor, SF_lep_sel[0], SF_lep_sel[1],SF_lep_sel[2], SF_lep_sel[3], err_lep_sel_up[0] if abs(tree.LepLepId[0])==11 else 0., err_lep_sel_up[1] if abs(tree.LepLepId[1])==11 else 0., err_lep_sel_up[2] if abs(tree.LepLepId[2])==11 else 0., err_lep_sel_up[3] if abs(tree.LepLepId[3])==11 else 0.))
                uncor_e_dn[ichan][k] += TRIG*(sigma_event(corr_factor, SF_lep_trig[0], SF_lep_trig[1], SF_lep_trig[2], SF_lep_trig[3], err_lep_trig_dn[0] if abs(tree.LepLepId[0])==11 else 0., err_lep_trig_dn[1] if abs(tree.LepLepId[1])==11 else 0., err_lep_trig_dn[2] if abs(tree.LepLepId[2])==11 else 0., err_lep_trig_dn[3] if abs(tree.LepLepId[3])==11 else 0.)) + SEL*(sigma_event(corr_factor, SF_lep_sel[0], SF_lep_sel[1],SF_lep_sel[2], SF_lep_sel[3], err_lep_sel_dn[0] if abs(tree.LepLepId[0])==11 else 0., err_lep_sel_dn[1] if abs(tree.LepLepId[1])==11 else 0., err_lep_sel_dn[2] if abs(tree.LepLepId[2])==11 else 0., err_lep_sel_dn[3] if abs(tree.LepLepId[3])==11 else 0.))
                uncor_mu_up[ichan][k] += TRIG*(sigma_event(corr_factor, SF_lep_trig[0], SF_lep_trig[1], SF_lep_trig[2], SF_lep_trig[3], err_lep_trig_up[0] if abs(tree.LepLepId[0])==13 else 0., err_lep_trig_up[1] if abs(tree.LepLepId[1])==13 else 0., err_lep_trig_up[2] if abs(tree.LepLepId[2])==13 else 0., err_lep_trig_up[3] if abs(tree.LepLepId[3])==13 else 0.)) + SEL*(sigma_event(corr_factor, SF_lep_sel[0], SF_lep_sel[1],SF_lep_sel[2], SF_lep_sel[3], err_lep_sel_up[0] if abs(tree.LepLepId[0])==13 else 0., err_lep_sel_up[1] if abs(tree.LepLepId[1])==13 else 0., err_lep_sel_up[2] if abs(tree.LepLepId[2])==13 else 0., err_lep_sel_up[3] if abs(tree.LepLepId[3])==13 else 0.))
                uncor_mu_dn[ichan][k] += TRIG*(sigma_event(corr_factor, SF_lep_trig[0], SF_lep_trig[1], SF_lep_trig[2], SF_lep_trig[3], err_lep_trig_dn[0] if abs(tree.LepLepId[0])==13 else 0., err_lep_trig_dn[1] if abs(tree.LepLepId[1])==13 else 0., err_lep_trig_dn[2] if abs(tree.LepLepId[2])==13 else 0., err_lep_trig_dn[3] if abs(tree.LepLepId[3])==13 else 0.)) + SEL*(sigma_event(corr_factor, SF_lep_sel[0], SF_lep_sel[1],SF_lep_sel[2], SF_lep_sel[3], err_lep_sel_dn[0] if abs(tree.LepLepId[0])==13 else 0., err_lep_sel_dn[1] if abs(tree.LepLepId[1])==13 else 0., err_lep_sel_dn[2] if abs(tree.LepLepId[2])==13 else 0., err_lep_sel_dn[3] if abs(tree.LepLepId[3])==13 else 0.))

            yield_e_up[ichan][k] += weight_nom/SF_tot_nom * SF_var_e_up
            yield_e_dn[ichan][k] += weight_nom/SF_tot_nom * SF_var_e_dn
            yield_mu_up[ichan][k] += weight_nom/SF_tot_nom * SF_var_mu_up
            yield_mu_dn[ichan][k] += weight_nom/SF_tot_nom * SF_var_mu_dn


    comb_e_up = []
    comb_e_dn = []
    comb_uncor_e_up = []
    comb_uncor_e_dn = []
    comb_mu_up = []
    comb_mu_dn = []
    comb_uncor_mu_up = []
    comb_uncor_mu_dn = []
    for ichan in range(len(chans)):
        comb_e_up.append(0.)
        comb_e_dn.append(0.)
        comb_uncor_e_up.append(0.)
        comb_uncor_e_dn.append(0.)
        comb_mu_up.append(0.)
        comb_mu_dn.append(0.)
        comb_uncor_mu_up.append(0.)
        comb_uncor_mu_dn.append(0.)

    # Sum relative uncertainties in quadrature
    for k in range (2):
        if ( correlated_leptons ):
            print "=================================================================="
            print "| CORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type)
            print "=================================================================="
            for ichan in range(len(chans)):
                comb_e_up[ichan] += ((yield_e_up[ichan][k]-nom_yield[ichan])/nom_yield[ichan]*100)**2
                comb_e_dn[ichan] += ((-yield_e_dn[ichan][k]+nom_yield[ichan])/nom_yield[ichan]*100)**2
                comb_mu_up[ichan] += ((yield_mu_up[ichan][k]-nom_yield[ichan])/nom_yield[ichan]*100)**2
                comb_mu_dn[ichan] += ((-yield_mu_dn[ichan][k]+nom_yield[ichan])/nom_yield[ichan]*100)**2
                print "{0}/ele = +{1:.1f}% -{2:.1f}%  {0}/mu = +{3:.1f}% -{4:.1f}%".format(chans[ichan], (yield_e_up[ichan][k]-nom_yield[ichan])/nom_yield[ichan]*100, (-yield_e_dn[ichan][k]+nom_yield[ichan])/nom_yield[ichan]*100, (yield_mu_up[ichan][k]-nom_yield[ichan])/nom_yield[ichan]*100, (-yield_mu_dn[ichan][k]+nom_yield[ichan])/nom_yield[ichan]*100)
            print "=================================================================="

        if (uncorrelated_leptons):
            print "=================================================================="
            print "| UNCORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type)
            print "=================================================================="
            for ichan in range(len(chans)):
                comb_uncor_e_up[ichan] += (uncor_e_up[ichan][k]/br_fs[ichan])
                comb_uncor_e_dn[ichan] += (uncor_e_dn[ichan][k]/br_fs[ichan])
                comb_uncor_mu_up[ichan] += (uncor_mu_up[ichan][k]/br_fs[ichan])
                comb_uncor_mu_dn[ichan] += (uncor_mu_dn[ichan][k]/br_fs[ichan])
                print "{0}/ele = +{1:.1f}% -{2:.1f}%  {0}/mu = +{3:.1f}% -{4:.1f}%".format(chans[ichan], 100*(uncor_e_up[ichan][k]/br_fs[ichan])**0.5,100*(uncor_e_dn[ichan][k]/br_fs[ichan])**0.5,100*(uncor_mu_up[ichan][k]/br_fs[ichan])**0.5,100*(uncor_mu_dn[ichan][k]/br_fs[ichan])**0.5)
            print "=================================================================="

    if ( correlated_leptons ):
        print "=========================================="
        print "| CORRELATED COMBINATION IN {} SAMPLE |".format(Type)
        print "=========================================="
        for ichan in range(len(chans)):
            print "{0}/ele = +{1:.1f}% -{2:.1f}%  {0}/mu = +{3:.1f}% -{4:.1f}%".format(chans[ichan], comb_e_up[ichan]**0.5,comb_e_dn[ichan]**0.5,comb_mu_up[ichan]**0.5,comb_mu_dn[ichan]**0.5)
        print "=============================================================================="

    if (uncorrelated_leptons):
        print "============================================"
        print "| UNCORRELATED COMBINATION IN {} SAMPLE |".format(Type)
        print "============================================"
        for ichan in range(len(chans)):
            print "{0}/ele = +{1:.1f}% -{2:.1f}%  {0}/mu = +{3:.1f}% -{4:.1f}%".format(chans[ichan], 100*(comb_uncor_e_up[ichan]**0.5),100*(comb_uncor_e_dn[ichan]**0.5),100*(comb_uncor_mu_up[ichan]**0.5),100*(comb_uncor_mu_dn[ichan]**0.5))
        print "=============================================================================="





