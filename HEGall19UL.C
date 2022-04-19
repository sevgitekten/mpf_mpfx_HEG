#include <iostream>
#define HEGall19UL_cxx
#include "HEGall19UL.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include <set>
#include "TRandom3.h"
#include <vector>
#include "basic.h"
#include <string>
#include <map>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include <time.h>

using namespace std;

double DPhi(double phi1, double phi2){
    double dphi_ = fabs(phi1-phi2);
    return (dphi_ <= TMath::Pi())? dphi_ : TMath::TwoPi() -dphi_;
}
double oplus(double a, double b) {
  return sqrt(a*a + b*b);
}

void HEGall19UL::Loop()
{
    time_t start,end; double dif;  time (&start);
    
    //===============================JSON control=======================================
     std::vector<std::vector<int>> vLS{
         //{67, 88},                               //306926  !! bu kesin kalkti
         {15, 213},                                 //306929  0
         {110, 1448},{1450, 1600},                  //306936  1 2
         {100, 1200},                               //307014  3
         {1, 397},                                  //307015  4
         {1, 317},                                  //307016  5
         {1, 761},                                  //307017  6
         {71, 573}, {575, 575}, {577, 577},         //307042  7 8 9
         {1, 86},                                   //307044  10
         {1, 17}, {19, 50},                         //307045  11 12
         {1, 526},                                  //307047  13
         {1, 191},                                  //307048  14
         {1, 24},                                   //307049  15
         {1, 320},                                  //307050  16
         {1, 513},                                  //307051  17
         {1, 24},                                   //307052  18
         {1, 57},                                   //307053  19
         {1, 62},                                   //307054  20
         {1, 878},                                  //307055  21
         {1, 795},                                  //307062  22
         {1, 230},                                  //307063  23
         {71, 388}, {845, 938},                     //307073  24 25
         {1, 620}, {622, 1669},                     //307076  26 27
         {110, 1300}                                //307082  28
         };
     
     /* Multimap with duplicates */
     std::multimap<int, int> mLS1 = {
        //{306926, 0 },
         {306929, 0 },
         {306936, 1 },{306936, 2},
         {307014, 3 },
         {307015, 4 },
         {307016, 5 },
         {307017, 6 },
         {307042, 7 },{307042, 8},{307042, 9},
         {307044, 10},
         {307045, 11},{307045, 12},
         {307047, 13},
         {307048, 14},
         {307049, 15},
         {307050, 16},
         {307051, 17},
         {307052, 18},
         {307053, 19},
         {307054, 20},
         {307055, 21},
         {307062, 22},
         {307063, 23},
         {307073, 24},{307073, 25},
         {307076, 26},{307076, 27},
         {307082, 28},
     };
     std::multimap<int, int> mLS2 (mLS1.begin(), mLS1.end());
      //===============================JSON control=======================================
     
    
    if (fChain == 0) return; //MakeClass
    Long64_t nentries = fChain->GetEntriesFast();//MakeClass
    Long64_t nbytes = 0, nb = 0; //MakeClass
    
    //cout<<"nentries"<<nentries<<endl;
    //nentries = 1000; //72372097 olay var!
    //TH2D *hotzonemap = (TH2D*)hotzone->Get("h2hotfilter");
    //TH2D *hotzonemap = (TH2D*)hotzone->Get("h2hot_ul17_plus_hbpw89");
    //TH2D *coldzonemap = (TH2D*)coldzone->Get("all/h2hole");
    TH2D *hot_coldzonemap = (TH2D*)hot_cold->Get("h2hot_hot_cold");
    TFile myFile("12Nisan_HEG_VX_mpf_133_6389_alleta.root", "RECREATE");
    
    
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V5_DATA_L1FastJet_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V5_DATA_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V5_DATA_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunH_VX_DATA_L2L3Residual_AK4PFchs.txt");
       
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    vParam_pfchs.push_back(*pfchs_l2l3res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    FactorizedJetCorrector *mpfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
    static const int netabins = 9;
       //static const double etabins[netabins+1] = {0,0.5};
       static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7,1.3};
       /*//static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
       //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
       //static const double etabins[netabins+1] = {0,1.3,3.2,4.7};
       //static const double etabins[netabins+1] = {0,1.3};
           
       };*/
    const double x[9][129]=
        {
            //=========================yeni Recobinninglerrrr=============================
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684, 1733, 1784, 1836, 1890, 1944, 2000, 2057, 2116, 2176, 2238, 2301, 2366, 2432, 2500, 2569, 2640, 2712, 2787, 2863, 2941, 3021, 3103, 3187, 3273, 3360, 3450, 3637, 3832, 4836, 6076, 6231, 6389}, // Eta_0.0-0.5
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684, 1733, 1784, 1836, 1890, 1944, 2000, 2057, 2116, 2176, 2238, 2301, 2366, 2432, 2500, 2569, 2640, 2712, 2787, 2863, 2941, 3021, 3103, 3187, 3273, 3450, 3637, 4364, 5220, 5355, 5492}, // Eta_0.5-1.0
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684, 1733, 1784, 1836, 1890, 1944, 2000, 2057, 2116, 2176, 2238, 2301, 2366, 2432, 2500, 2569, 2640, 2787, 2941, 3360, 3832}, // Eta_1.0-1.5
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684, 1733, 1784, 1836, 1890, 1944, 2000, 2057, 2116, 2301, 2500, 2569, 2640}, // Eta_1.5-2.0
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684}, // Eta_2.0-2.5
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032}, // Eta_2.5-3.0
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032}, // Eta_3.0-3.2
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032}, // Eta_3.2-4.7
            {10, 11, 12, 13.5, 15, 16.5, 18, 19.5, 21, 22.5, 24, 26, 28, 30, 32, 34.4, 37, 39.9, 43, 46.1, 49, 52.5, 56, 60, 64, 69, 74, 79, 84, 90, 97, 105, 114, 123, 133, 143, 153, 163, 174, 185, 196, 208, 220, 232, 245, 258, 272, 286, 300, 315, 330, 346, 362, 378, 395, 412, 430, 449, 468, 487, 507, 527, 548, 570, 592, 614, 638, 662, 686, 711, 737, 763, 790, 818, 846, 875, 905, 935, 967, 999, 1032, 1066, 1101, 1136, 1172, 1209, 1248, 1287, 1327, 1368, 1410, 1453, 1497, 1542, 1588, 1636, 1684, 1733, 1784, 1836, 1890, 1944, 2000, 2057, 2116, 2176, 2238, 2301, 2366, 2432, 2500, 2569, 2640, 2712, 2787, 2863, 2941, 3021, 3103, 3187, 3273, 3360, 3450, 3637, 3832, 4836, 6076, 6231, 6389}, // Eta_0.0-1.3
        };
      
       const int nx[9] = {128,126,116,108,96,80,80,80,128};
    
    // const int nx[1] = {64};
    // const int nx[3] = {64,63,58};
    // const int nx[8] = {64,63,58,54,48,40,40,40};
    //const int nx[9] = {64,63,58,54,48,40,40,40,64};
    //const int nx[9] = {127,126,116,108,96,80,80,80,127};
    //const int nx[1] = {128};
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    //std::vector<TH1F*> vectptt[7];
    
     /*std::vector<TProfile*> vectptt; //dataquality icin
     std::vector< std::vector< TProfile*> > V2daq;   //dataquality icin
     std::vector< std::vector< std::vector< TProfile*> > > VVdaq;    //dataquality icin
    */
    std::vector<TH2D*> vectptt; //mpf icin
    std::vector< std::vector< TH2D*> > V2daq;   //mpf icin
    std::vector< std::vector< std::vector< TH2D*> > > VVdaq;
    
    const int n3 = 500;
    vector<double> v3(n3+1);
    for (unsigned int i = 0; i != n3+1; ++i) v3[i] = -5. + i*0.02;
    
    for(int l=0;l<6;++l){ //triggerlara göre, dataqualityde 2 boyutlu vector cizdirdiğinde  acarsin
    for ( int k = 0; k<netabins; ++k)
    {
       /*TProfile *pchf_0= new TProfile(Form("pchf_%d_%d",l,k),"pchf",nx[k],&x[k][0]); vectptt.push_back (pchf_0); pchf_0->Sumw2();
       TProfile *pnhf_1= new TProfile(Form("pnhf_%d_%d",l,k),"pnhf",nx[k],&x[k][0]); vectptt.push_back (pnhf_1); pnhf_1->Sumw2();
       TProfile *pnef_2= new TProfile(Form("pnef_%d_%d",l,k),"pnef",nx[k],&x[k][0]); vectptt.push_back (pnef_2); pnef_2->Sumw2();
       TProfile *pcef_3= new TProfile(Form("pcef_%d_%d",l,k),"pcef",nx[k],&x[k][0]); vectptt.push_back (pcef_3); pcef_3->Sumw2();
       TProfile *pmuf_4= new TProfile(Form("pmuf_%d_%d",l,k),"pmuf",nx[k],&x[k][0]); vectptt.push_back (pmuf_4); pmuf_4->Sumw2();
       TProfile *ppuf_5= new TProfile(Form("ppuf_%d_%d",l,k),"ppuf",nx[k],&x[k][0]); vectptt.push_back (ppuf_5); ppuf_5->Sumw2();
       TProfile *phhf_6= new TProfile(Form("phhf_%d_%d",l,k),"phhf",nx[k],&x[k][0]); vectptt.push_back (phhf_6); phhf_6->Sumw2();
       TProfile *phef_7= new TProfile(Form("phef_%d_%d",l,k),"phef",nx[k],&x[k][0]); vectptt.push_back (phef_7); phef_7->Sumw2();
       
       
       TProfile *pchftp_8  = new TProfile(Form("pchftp_%d_%d",l,k),"pchftp",nx[k],&x[k][0]); vectptt.push_back (pchftp_8);  pchftp_8->Sumw2();
       TProfile *pnhftp_9  = new TProfile(Form("pnhftp_%d_%d",l,k),"pnhftp",nx[k],&x[k][0]); vectptt.push_back (pnhftp_9);  pnhftp_9->Sumw2();
       TProfile *pneftp_10 = new TProfile(Form("pneftp_%d_%d",l,k),"pneftp",nx[k],&x[k][0]); vectptt.push_back (pneftp_10); pneftp_10->Sumw2();
       TProfile *pceftp_11 = new TProfile(Form("pceftp_%d_%d",l,k),"pceftp",nx[k],&x[k][0]); vectptt.push_back (pceftp_11); pceftp_11->Sumw2();
       TProfile *pmuftp_12 = new TProfile(Form("pmuftp_%d_%d",l,k),"pmuftp",nx[k],&x[k][0]); vectptt.push_back (pmuftp_12); pmuftp_12->Sumw2();
       TProfile *ppuftp_13 = new TProfile(Form("ppuftp_%d_%d",l,k),"ppuftp",nx[k],&x[k][0]); vectptt.push_back (ppuftp_13); ppuftp_13->Sumw2();
       TProfile *phhftp_14 = new TProfile(Form("phhftp_%d_%d",l,k),"phhftp",nx[k],&x[k][0]); vectptt.push_back (phhftp_14); phhftp_14->Sumw2();
       TProfile *pheftp_15 = new TProfile(Form("pheftp_%d_%d",l,k),"pheftp",nx[k],&x[k][0]); vectptt.push_back (pheftp_15); pheftp_15->Sumw2();
        */
        TH2D *h2mpf_0= new TH2D (Form("h2mpf_%d_%d",l,k),"h2mpf",nx[k],&x[k][0],n3,&v3[0]); vectptt.push_back (h2mpf_0); h2mpf_0->Sumw2();
        TH2D *h2mpfx_1= new TH2D (Form("h2mpfx_%d_%d",l,k),"h2mpfx",nx[k],&x[k][0],n3,&v3[0]); vectptt.push_back (h2mpfx_1); h2mpfx_1->Sumw2();
        
       
        V2daq.push_back(vectptt);vectptt.clear();
    }
    VVdaq.push_back(V2daq);V2daq.clear(); //2 boyutlu vektoru olusturdugunda ac
    } //l nin yani triggerların
    TLorentzVector p4pf;
    TLorentzVector mp4pf;
    
    //trigger PT aralıklarını buraya yerlestırdım
   // vector<double>vecttrigcuts={43,64,64,84,114,272,6389}; //15Aralik
    vector<double>vecttrigcuts={74,74,74,97,133,272,6389};   //22 Kasim yeni trigger rangeler - 19Aralik
    //vector<double>vecttrigcuts={49,49,64,97,133,272,6389};   //yeni trigger rangeler
    
   
    //=================Define the variables for Correction========================
    
    float pf_jtpt[100],pf_jt_eta[100],pf_jt_phi[100];    //corr definition (disarda kalabilir)
    vector<float>pfchs_ptcorr;
    
    float mpf_jtpt[100],mpf_jt_eta[100],mpf_jt_phi[100];    //corr definition (disarda kalabilir)
    vector<float>mpfchs_ptcorr;
    
    //---------------------------events dongusu icine girdik---------
    
    for (Long64_t jentry=0; jentry<nentries;jentry++)   //MakeClass
    {
        Long64_t ientry = LoadTree(jentry);             //MakeClass
        if (ientry < 0) break;                          //MakeClass
        nb = fChain->GetEntry(jentry);   nbytes += nb;  //MakeClass
        
       // cout<<"\r"<<"event number: "<<jentry<<"/ "<<nentries<<flush;
        
        if (TriggerDecision_.empty()) continue;
        //double l_pt1=0.0, l_pt2=0.0, l_pt3=0.0;
        if (FilterDecision_.size()!=0)  continue;
        
        int pf_njt =PFJetsCHS__;        //==========corr definition
        pfchs_ptcorr.resize(pf_njt);    //=======corr definition
        
        // tag and probe  jets variables
         double jet0_pt=0; double jet1_pt=0; double jet2_pt=0;
         double jet0_phi=0; double jet1_phi=0; double jet2_phi=0;
         double jet0_eta=99; double jet1_eta=99; double jet2_eta=99;
         //double dphi=0; double ptave=0; double alpha=1;

         
         //find leading jets (JEC may change ordering)
         //less then 3 jets we input -1 on the place of the index
         Int_t           jt3leads[3]; // The three leading jets
         Double_t        jtpt[PFJetsCHS__];
         Double_t        mjtpt[PFJetsCHS__];
         Double_t        jteta[PFJetsCHS__];
         Double_t        jtphi[PFJetsCHS__];
         for (int i = 0; i<3; ++i) jt3leads[i] = -1;
         vector<int> qcdpass4;
         vector<int> qcdpass5;
        
        // Variables for MPF
        double metsumet = PFMet__sumEt_;
        double met    = PFMet__et_;
        double metphi = PFMet__phi_;
        double mex    = met*cos(metphi);
        double mey    = met*sin(metphi);
        double metsumet1 = metsumet;
         
        
     
       Int_t           mjt3leads[3]; // The three leading jets
       for (int i = 0; i<3; ++i) mjt3leads[i] = -1;
       /////////////////////////////////////////////////////
       ///Jet loop for  mpf?
       /////////////////////////////////////////////////////
        
         
         for (int j=(PFJetsCHS__-1); j>-1; j--)
              {
                  
                  mp4pf.SetPxPyPzE(PFJetsCHS__P4__fCoordinates_fX[j],PFJetsCHS__P4__fCoordinates_fY[j],
                                  PFJetsCHS__P4__fCoordinates_fZ[j],PFJetsCHS__P4__fCoordinates_fT[j]);
         
              //================Calculation of correction ///
                 // cout<<"mpftest0 "<<endl;
                  mpf_jtpt[j]   = mp4pf.Pt(); mpf_jt_eta[j] = mp4pf.Eta(); mpf_jt_phi[j] = mp4pf.Phi();
                  //double pf_pt_old=mpf_jtpt[j];
                  
                  float rho=EvtHdr__mPFRho;
                  double area=PFJetsCHS__area_[j];
                  double ptraw = mpf_jtpt[j]/PFJetsCHS__cor_[j]; //undo islemi
                 
                  mpfchs_jec->setJetEta(mpf_jt_eta[j]);mpfchs_jec->setJetPt(ptraw);
                  mpfchs_jec->setRho(rho);mpfchs_jec->setJetA(area);
                  
                  mpfchs_ptcorr[j]= ptraw*mpfchs_jec->getCorrection();
                  
                  double pf_pt=mpfchs_ptcorr[j]; double pf_jteta=mpf_jt_eta[j];double pf_phi=mpf_jt_phi[j];
                  mjtpt[j] = pf_pt; //jteta[j]= pf_jteta; jtphi[j]= pf_phi; //tag and probe corr
                  
              //=====================the end of correction////////

         //cout<<"mpftest01 "<<endl;
                  //==============finding leading jets
                    if (mjt3leads[0]==-1 or mjtpt[mjt3leads[0]]<mjtpt[j]) {
                      mjt3leads[2] = mjt3leads[1];
                      mjt3leads[1] = mjt3leads[0];
                      mjt3leads[0] = j;
                    } else if (mjt3leads[1]==-1 or mjtpt[mjt3leads[1]]<mjtpt[j]) {
                      mjt3leads[2] = mjt3leads[1];
                      mjt3leads[1] = j;
                    } else if (mjt3leads[2]==-1 or mjtpt[mjt3leads[2]]<mjtpt[j]) {
                      mjt3leads[2] = j;
                    }
                 
                  
              }
         //cout<<"mpftest02 "<<endl;
         // qcdpass2: 1 2  qcdpass3: 0 2
          bool mpass1=false;
          bool mpass2=false;
          bool mpass3=false;
              if(mjt3leads[0]>-1) mpass1=true;
              if(mjt3leads[1]>-1) mpass2=true;
              if(mjt3leads[2]>-1) mpass3=true;
              //if (qcdpass.size()>=2) cout<<" jt3leads[0]  :"<<jt3leads[0] <<" jt3leads[1]  :"<<jt3leads[1]<<" jt3leads[2]  :"<<jt3leads[2]<<" qcdpass: "<<qcdpass[jj]<<"   jentry: "<<jentry<<endl;

         // cout<<"mpftest02 "<<endl;
          
          //-----------------------------------Tag and probe-----------
          //the leading indices
          
          int mi0 = mjt3leads[0]; int mi1 = mjt3leads[1]; int mi2 = mjt3leads[2];
          if (mi0 < 0. ) return; //This should not happen
         double mptave=0.;
          if(mpass1 and mpass2){
              mptave = (mi1>=0 ? 0.5 * (mjtpt[mi0] + mjtpt[mi1]) : mjtpt[mi0]);
          }
        // cout<<"mpftest1 "<<endl;

        
        
        
         /////////////////////////////////////////////////////
         ///Jet Loop
         ////////////////////////////////////////////////////
       
        for (int j=(PFJetsCHS__-1); j>-1; j--)
        //for(int j=0; j<PFJetsCHS__; j++)
        {
            
            p4pf.SetPxPyPzE(PFJetsCHS__P4__fCoordinates_fX[j],PFJetsCHS__P4__fCoordinates_fY[j],
                            PFJetsCHS__P4__fCoordinates_fZ[j],PFJetsCHS__P4__fCoordinates_fT[j]);
            
            
            ////////////correction icin eklediklerimmm///
            pf_jtpt[j]   = p4pf.Pt(); pf_jt_eta[j] = p4pf.Eta(); pf_jt_phi[j] = p4pf.Phi();
            double pf_pt_old=pf_jtpt[j];
            
            float rho=EvtHdr__mPFRho;
            double area=PFJetsCHS__area_[j];
            double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j];
            //double ptraw = pf_jtpt[j]+rho*area; //undo yapmadan ptraw hesabi
            pfchs_jec->setJetEta(pf_jt_eta[j]);
            pfchs_jec->setJetPt(ptraw);
            pfchs_jec->setRho(rho);
            pfchs_jec->setJetA(area);
            
            pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
            double pf_pt=pfchs_ptcorr[j];
            double pf_jteta=pf_jt_eta[j];
            double pf_phi=pf_jt_phi[j];
            jtpt[j] = pf_pt; jteta[j]= pf_jteta; jtphi[j]= pf_phi; //tag and probe corr
                
            //=====================the end of correction////////
            
                //====================mpf type1met=============
                //only use jets with corrected pt>Recopt to equalize DATA and MC thresholds
                if(fabs(pf_jteta)<4.7){
                    if(pf_pt_old > 15.){
                        double dpt= -pf_pt +pf_pt_old;
                        mex += dpt * cos(pf_phi);
                        mey += dpt * sin(pf_phi);
                        metsumet1 += pf_pt -pf_pt_old;
                    }
                }
            
                //==============finding leading jets
                  if (jt3leads[0]==-1 or jtpt[jt3leads[0]]<jtpt[j]) {
                    jt3leads[2] = jt3leads[1];
                    jt3leads[1] = jt3leads[0];
                    jt3leads[0] = j;
                  } else if (jt3leads[1]==-1 or jtpt[jt3leads[1]]<jtpt[j]) {
                    jt3leads[2] = jt3leads[1];
                    jt3leads[1] = j;
                  } else if (jt3leads[2]==-1 or jtpt[jt3leads[2]]<jtpt[j]) {
                    jt3leads[2] = j;
                  }
                   
            
            if( PFJetsCHS__tightID_[j] && (PFMet__et_<0.3*PFMet__sumEt_) && (PFJetsCHS__cemf_[j]<0.9) && (PFJetsCHS__muf_[j]<0.9))
            {
                char name[100];
                sprintf(name,"%d",EvtHdr__mRun);
                for ( int k = 0; k<netabins; ++k)
                {
                    double etamin=etabins[k];double etamax=etabins[k+1];
                    if (k==8){etamin=etabins[0]; etamax=etabins[9];}
                    
                    if (fabs(pf_jteta)>= etamin && fabs(pf_jteta)< etamax)
                    //if ( fabs(pf_jteta)>=etabins[k] && fabs(pf_jteta)<etabins[k+1])
                    {
                        if(mptave>0.){
                        for(int trnameindex=4; trnameindex<6 ; trnameindex++)
    
                        {
                           
                            for(int trdecindex=0; trdecindex<TriggerDecision_.size(); trdecindex++)
                            {
                                 
                               // if(TriggerDecision_[trdecindex]==trnameindex && ((pf_pt)>=vecttrigcuts[trnameindex]&& (pf_pt)<vecttrigcuts[trnameindex+1]))
                                if(TriggerDecision_[trdecindex]==trnameindex && ((mptave)>=vecttrigcuts[trnameindex] && (mptave)<vecttrigcuts[trnameindex+1]))
                                {
                                                            if (mLS1.find(EvtHdr__mRun)!=mLS1.end()){ //======
                                                            //cout<<"Run Number in the JSON " <<EvtHdr__mRun<<endl;
                                                            int loop_control=0; //======
                                                            auto LS_range = mLS1.equal_range(EvtHdr__mRun);//======
                                                            for (auto it=LS_range.first; it != LS_range.second; ++it)//======
                                                                {//======
                                                                    
                                                                    if (EvtHdr__mLumi>=vLS[it->second][0] && EvtHdr__mLumi<=vLS[it->second][1] )//======
                                                                    { loop_control=loop_control+1;//======
                                                            //cout <<" in the JSON " << trnameindex<< "  "<<EvtHdr__mRun <<"  "<< EvtHdr__mLumi <<endl; //======
                                                            int _prescale =1;
                                                            int l1_ps= L1Prescale_[trdecindex]; //L1 Trigger
                                                            int hlt_ps= HLTPrescale_[trdecindex]; //HLT Trigger
                                                            
                                                            if (l1_ps == 0) l1_ps =1;
                                                            if (hlt_ps == 0) hlt_ps =1;
                                                            _prescale = l1_ps * hlt_ps;
                                                            if (_prescale>=1 && _prescale<1000000){
                                                                
                                                                //int hotzonebin=hotzonemap->FindBin(pf_jteta,pf_phi);
                                                                //int coldzonebin=coldzonemap->FindBin(pf_jteta,pf_phi);
                                                                int hot_coldzonebin=hot_coldzonemap->FindBin(pf_jteta,pf_phi);
                                                                //if(hotzonemap->GetBinContent(hotzonebin)==0.0 && (coldzonemap->GetBinContent(coldzonebin)==0.0))
                                                                if(hot_coldzonemap->GetBinContent(hot_coldzonebin)==0.0)
                                                                {
                                                                    
                                                                /*VVdaq[trnameindex][k][0]->Fill(pf_pt,PFJetsCHS__chf_[j]);
                                                                VVdaq[trnameindex][k][1]->Fill(pf_pt,PFJetsCHS__nhf_[j]);
                                                                VVdaq[trnameindex][k][2]->Fill(pf_pt,PFJetsCHS__nemf_[j]);
                                                                VVdaq[trnameindex][k][3]->Fill(pf_pt,PFJetsCHS__cemf_[j]);
                                                                VVdaq[trnameindex][k][4]->Fill(pf_pt,PFJetsCHS__muf_[j]);
                                                                VVdaq[trnameindex][k][5]->Fill(pf_pt,PFJetsCHS__betaPrime_[j]);
                                                                VVdaq[trnameindex][k][6]->Fill(pf_pt,PFJetsCHS__hf_hf_[j]);
                                                                VVdaq[trnameindex][k][7]->Fill(pf_pt,PFJetsCHS__hf_phf_[j]);
                                                                   */
                                                                    if((k==8) && (trnameindex==4) ) qcdpass4.push_back(j);
                                                                    if((k==8) && (trnameindex==5) ) qcdpass5.push_back(j);
                                                                                                                }//hotzone
                                                                                                        }  //prescale>0
                                                                                                } //lumisection ls cut from JSON //======
                                                                        
                                                                                    } //Run number iteretive loop for JSON //======
                                                                        //if(loop_control==0){cout <<"out of the JSON " << EvtHdr__mRun <<" , "<< EvtHdr__mLumi <<endl;} //======
                                                                                } //Run check //======
                                                           /* else { //======
                                                                cout<<"out of the JSON " <<EvtHdr__mRun<<endl; //======
                                                            } //====== */
                                } //TriggerDecision
                            } //trdecindex
                        } //trnameindex
                      }//mptave
                    } // fabs
                } //etabin
            } //tightid
        } //jets
        
        double met1   = oplus(mex,mey);
        double metphi1= atan2(mey,mex);
        
        bool pass1=false;
        bool pass2=false;
        bool pass3=false;
        for(unsigned int jj=0; jj<qcdpass4.size(); ++jj){
            if(qcdpass4[jj]==jt3leads[0]) pass1=true;
            if(qcdpass4[jj]==jt3leads[1]) pass2=true;
            if(qcdpass4[jj]==jt3leads[2]) pass3=true;
            //if (qcdpass.size()>=2) cout<<" jt3leads[0]  :"<<jt3leads[0] <<" jt3leads[1]  :"<<jt3leads[1]<<" jt3leads[2]  :"<<jt3leads[2]<<" qcdpass: "<<qcdpass[jj]<<"   jentry: "<<jentry<<endl;
        }
        for(unsigned int jj=0; jj<qcdpass5.size(); ++jj){
            if(qcdpass5[jj]==jt3leads[0]) pass1=true;
            if(qcdpass5[jj]==jt3leads[1]) pass2=true;
            if(qcdpass5[jj]==jt3leads[2]) pass3=true;
        }
        
        
        //-----------------------------------Tag and probe-----------
        //the leading indices
        
        int i0 = jt3leads[0]; int i1 = jt3leads[1]; int i2 = jt3leads[2];
        if (i0 < 0. ) return; //This should not happen
        
        if(pass1 and pass2){
            
        double ptave = (i1>=0 ? 0.5 * (jtpt[i0] + jtpt[i1]) : jtpt[i0]);
        double dphi = (i1>=0 ? DPhi(jtphi[i0], jtphi[i1]) : 0.);
        double dpt = (i1>=0 ? fabs(jtpt[i0]-jtpt[i1])/(2*ptave) : 0.999);
        // If the jetID is bad for the third jet (and the third jet is visible), we set pt3 to ptave (alpha = 1)
        double pt3 = ((i1>=0 and i2>=0 and pass3 and jtpt[i2]>15) ? ((PFJetsCHS__>=3) ? jtpt[i2] : ptave) : 0.);
        double alpha = pt3/ptave;
        //cout<<endl<<"evnt 1  "<<EvtHdr__mEvent <<"  ,mptave  "<<mptave<< "  ,mjtpt[mi0],  "<<mjtpt[mi0]<< "  ,mjtpt[mi1],  "<<mjtpt[mi1]<< "  ,mjtpt[mi2],  "<<mjtpt[mi2]<<endl;

        //cout<<endl<<"evnt 2  "<<EvtHdr__mEvent <<"  ,ptave  "<<ptave<< "  ,jtpt[i0],  "<<jtpt[i0]<< "  ,jtpt[i1],  "<<jtpt[i1]<< "  ,jtpt[i2],  "<<jtpt[i2]<<endl;


        if ( i0>=0 and jtpt[i0]>15 and fabs(jteta[i0])<1.3 ) { // First leading jet
            if (i1>=0 and jtpt[i1]>15 and fabs(jteta[i1])<1.3 ) { // Second leading jet
               if (alpha < 1.0) { //first loose alpha cut for mpf
                  for (auto itag_lead = 0u; itag_lead<2u; ++itag_lead) { // Look for both t&p combos for the leading jets
                       int itag = jt3leads[itag_lead];
                       int iprobe = jt3leads[(itag_lead==0 ? 1 : 0)];
                       double etatag = jteta[itag];
                       double etaprobe = jteta[iprobe];
                       double pttag = jtpt[itag];
                       double ptprobe = jtpt[iprobe];
                       double phiprobe = jtphi[iprobe];
                        
                        // Dijet balance
                        if ((alpha < 0.3) && (dphi > 2.7)) { // Back-to-back condition
                            /*for(unsigned int t=0; t<qcdpass4.size(); ++t){
                                if(itag==qcdpass4[t]){
                                    VVdaq[4][8][8]->Fill(pttag,PFJetsCHS__chf_[iprobe]);
                                    VVdaq[4][8][9]->Fill(pttag,PFJetsCHS__nhf_[iprobe]);
                                    VVdaq[4][8][10]->Fill(pttag,PFJetsCHS__nemf_[iprobe]);
                                    VVdaq[4][8][11]->Fill(pttag,PFJetsCHS__cemf_[iprobe]);
                                    VVdaq[4][8][12]->Fill(pttag,PFJetsCHS__muf_[iprobe]);
                                    VVdaq[4][8][13]->Fill(pttag,PFJetsCHS__betaPrime_[iprobe]);
                                    VVdaq[4][8][14]->Fill(pttag,PFJetsCHS__hf_hf_[iprobe]);
                                    VVdaq[4][8][15]->Fill(pttag,PFJetsCHS__hf_phf_[iprobe]);
                                }
                            }

                            for(unsigned int t=0; t<qcdpass5.size(); ++t){
                                if(itag==qcdpass5[t]){
                                    VVdaq[5][8][8]->Fill(pttag,PFJetsCHS__chf_[iprobe]);
                                    VVdaq[5][8][9]->Fill(pttag,PFJetsCHS__nhf_[iprobe]);
                                    VVdaq[5][8][10]->Fill(pttag,PFJetsCHS__nemf_[iprobe]);
                                    VVdaq[5][8][11]->Fill(pttag,PFJetsCHS__cemf_[iprobe]);
                                    VVdaq[5][8][12]->Fill(pttag,PFJetsCHS__muf_[iprobe]);
                                    VVdaq[5][8][13]->Fill(pttag,PFJetsCHS__betaPrime_[iprobe]);
                                    VVdaq[5][8][14]->Fill(pttag,PFJetsCHS__hf_hf_[iprobe]);
                                    VVdaq[5][8][15]->Fill(pttag,PFJetsCHS__hf_phf_[iprobe]);
                                }
                            }
                          */
                        } //alpha ve dphi cut
                      
                      for(unsigned int t=0; t<qcdpass4.size(); ++t){
                          if(itag==qcdpass4[t] ){ //ptave cutini 11 Nisanda ekledim
                          //define mpf and mpfx
                          double mpf    = met1*cos(metphi1 - jtphi[itag])/2;
                          double mpfx   = met1*sin(metphi1 - jtphi[itag])/2;
                          mpf /= ptave;
                          mpfx /= ptave;
                          //fill mpf and mpfx
                          VVdaq[4][8][0]->Fill(ptave,mpf,1);
                          VVdaq[4][8][1]->Fill(ptave,mpfx,1);
                          //VVdaq[2][8][8]->Fill(pttag,PFJetsCHS__chf_[iprobe]);
                          //VVdaq[2][8][9]->Fill(pttag,PFJetsCHS__nhf_[iprobe]);
                            }
                      }
                      
                      for(unsigned int t=0; t<qcdpass5.size(); ++t){
                      if(itag==qcdpass5[t] ){ //ptave cutini 11 Nisanda ekledim
                          //define mpf and mpfx
                          double mpf    = met1*cos(metphi1 - jtphi[itag])/2;
                          double mpfx   = met1*sin(metphi1 - jtphi[itag])/2;
                          mpf /= ptave;
                          mpfx /= ptave;
                          //fill mpf and mpfx
                          VVdaq[5][8][0]->Fill(ptave,mpf,1);
                          VVdaq[5][8][1]->Fill(ptave,mpfx,1);
                          //VVdaq[2][8][8]->Fill(pttag,PFJetsCHS__chf_[iprobe]);
                          //VVdaq[2][8][9]->Fill(pttag,PFJetsCHS__nhf_[iprobe]);
                            }
                      }
                      

                      // cout<<"   itag    :"<<itag<<"     iprobe :    "<<iprobe<< "    Event   :" <<jentry <<"    etatag  :"<< etatag <<"   etaprobe  :   "<<etaprobe<<endl;
                      
                  }

               }
                
            }
             
         }
        }
        
    } //events
    
    myFile.cd();
    myFile.Write();
    myFile.Close();
    hotzone->Close();
    
    time (&end); dif = difftime (end,start); cout<< endl<< "zaman:"<< dif<<endl;
}


