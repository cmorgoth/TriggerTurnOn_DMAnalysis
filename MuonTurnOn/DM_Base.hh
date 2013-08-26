#ifndef DM_BASE_HH
#define DM_BASE_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <vector>
#include <map>
#include  <string>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include "TLorentzVector.h"

//Base Dark Matter Analysis Class
//Every other class (for each individual bkg channel and Data ) inherits from this one
//Defines de metods and constants commom for all the channels
//Baseline cuts, plots for MR and RSQ in different signal and control samples

class BaseDM{
  
public:
  
  static const int MR_Bins = 3;
  static const int RSQ_Bins = 4;
  
  static const float Lumi = 1.5*2.5*5.;//fb^-1
  
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];
  
  static const int btagIndex = 4;//0->Veto Btag(Loose), 1-> Btag(Loose) >=1, 2-> BtagTight (Loose, tight), 4(Med and Tight)
  
  BaseDM();
  BaseDM( const char*, TString , int);//constructor for Data
  BaseDM( const char*, TString, float, int );//constructor for MC
  
  ~BaseDM();
  
  virtual TH1F PlotMR_2Box();
  virtual TH1F PlotMR_1Box();
  virtual TH1F PlotMR_0Box();
  
  virtual TH1F PlotRSQ_2Box();
  virtual TH1F PlotRSQ_1Box();
  virtual TH1F PlotRSQ_0Box();
  
  virtual TH2F PlotRSQ_vs_MR_0Box( );
  virtual TH2F PlotRSQ_vs_MR_1Box( );
  virtual TH2F PlotRSQ_vs_MR_2Box( );

  virtual bool pfJetPassCSVM(double );
  virtual int pfJetPassCSVM(double*, int);
  virtual std::vector<TH2F*> Plot_2DRazor();
  
  virtual bool CalcWeight();
  virtual bool PrintEvents();
  virtual double GetWeight(){return weight;};

  virtual std::vector<TH1F*> PlotMETx();
  virtual std::vector<TH1F*> PlotMETy();
  virtual std::vector<TH1F*> PlotMETmag();
  virtual std::vector<TH1F*> PlotMETphi(){};
  virtual std::vector<TH1F*> DoubleMuBoxPlots();
  
  virtual std::vector<TH1F*> PlotHT();
  virtual std::vector<TH1F*> PlotMHTmag(){};
  virtual std::vector<TH1F*> PlotMHTphi(){};
  
  virtual bool SetBrachStatus();

  virtual double HLTscale(double, double);
  virtual double HLTscaleEle(double, double);

  // a:Loose Btag, b:Medium Btag, c: Tight Btag                                 
  virtual bool SetBtagCut(int a, int b, int c){nBtagCut[0]=a; nBtagCut[1]=b; nBtagCut[2]=c;};
  
private:
  TTree* T;
  TFile* F;
  
  TH2F* hlt;
  TH2F* hlt_ele;

  TEfficiency* eff;
  TEfficiency* eff_ele;
  
  int metIndex;//0, 1, 2, 3 (2: type1 pfMET Correction)
  float MRMin;
  float RSQMin;
  float weight;
  float sigma;
  TString processName;
  bool fBtag[5];
  int nBtagCut[3];
  TString BtagBranch;

};



#endif
