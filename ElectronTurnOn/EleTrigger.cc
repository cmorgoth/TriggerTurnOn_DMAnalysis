#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <iostream>
#include "DM_2DRatio.hh"
#include "DM_Base.hh"


//const float BaseDM::RSQ_BinArr[] = {0.3, 0.4, 0.5, 0.6, 0.8, 2.5};//Mu
//const float BaseDM::MR_BinArr[] = {200., 300., 400. ,500., 3500.};//Mu

const float BaseDM::RSQ_BinArr[] = {0.3, 0.4, 0.5, 0.6, 2.5};//Mu                                                   
const float BaseDM::MR_BinArr[] = {200., 300., 400. , 3500.};//Mu

int scaleFactor(double MR, double RSQ);
double FindBin(double MR, double R2);

//TFile* f;
TH2F* num1;
TEfficiency* pEff2;

int main(){

  gROOT->Reset();
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  //const double RSQ_BinArr[] = {0.3, 0.4, 0.5, 0.6, 0.8, 2.5};//Mu
  //const double MR_BinArr[] = {200., 300., 400. ,500., 3500.};//Mu
  const double RSQ_BinArr[] = {0.3, 0.4, 0.5, 0.6, 2.5};//Ele
  const double MR_BinArr[] = {200., 300., 400., 3500.};//Ele

  TFile* f = new TFile("/media/data/cmorgoth/Data/DMData/TriggerData/SingleEle_ILV_ABCD_NOELEVETO.root");
  TTree* t = (TTree*)f->Get("outTree");
  
  double MR[4], RSQ[4];
  int HLT_Razor_prescaled, HLT_Razor, BOX;
  
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("MR",1);
  t->SetBranchStatus("RSQ",1);
  t->SetBranchStatus("HLT_Razor_prescaled", 1);
  t->SetBranchStatus("HLT_Razor",1);
  t->SetBranchStatus("BOX_NUM",1);
  
  t->SetBranchAddress("MR", MR);
  t->SetBranchAddress("RSQ", RSQ);
  t->SetBranchAddress("HLT_Razor", &HLT_Razor);
  t->SetBranchAddress("HLT_Razor_prescaled", &HLT_Razor_prescaled);
  t->SetBranchAddress("BOX_NUM", &BOX);
  
  num1 = new TH2F("num1","num1", 3, MR_BinArr, 4, RSQ_BinArr);
  TH2F* den1 = new TH2F("den1","den1", 3, MR_BinArr, 4, RSQ_BinArr);
  
  TH1F* MRn1 = new TH1F("MRn1","Elehad-Turn-On-MR", 3, MR_BinArr);
  TH1F* MRd1 = new TH1F("MRd1","Elehad-Turn-On-MR-Den", 3, MR_BinArr);

  TH1F* R2n1 = new TH1F("R2n1","Elehad-Turn-On-R2", 4, RSQ_BinArr);
  TH1F* R2d1 = new TH1F("R2d1","Elehad-Turn-On-R2-Den", 4, RSQ_BinArr);

  //std::cout << "Nentries: " << t->GetEntries() << std::endl;
  
  for(int i = 0; i < t->GetEntries(); i++ ){
    t->GetEntry(i);
    if( RSQ[2] >= 0.30 && MR[2] >= 200. && BOX == 0){
      den1->Fill(MR[2], RSQ[2]);
      MRd1->Fill(MR[2]);
      R2d1->Fill(RSQ[2]);
      if( HLT_Razor == 1 ){
	num1->Fill(MR[2], RSQ[2]);
	MRn1->Fill(MR[2]);
	R2n1->Fill(RSQ[2]);
      }
    }
    
  }
  
  num1->Sumw2();
  den1->Sumw2();
  pEff2 = new TEfficiency(*num1, *den1);
  std::cout << "Num: " << num1->Integral() << std::endl;
  std::cout << "Den: " << den1->Integral() << std::endl;
  den1->SetStats(0);
  num1->Divide(den1);
  num1->SetNameTitle("h", "h1");
  num1->SetStats(0);
  num1->Draw("colz");
  C1->SaveAs("SingleElectronPD_2d_Trigger.png");
  C1->SaveAs("SingleElectronPD_2d_Trigger.pdf");

  pEff2->Draw("colz");
  C1->SaveAs("EleEff_2d.pdf");
  C1->SaveAs("EleEff_2d.png");

  TString s;
  TString s1;
  
  MRn1->Sumw2();
  MRd1->Sumw2();
  TEfficiency* pEff3 = new TEfficiency(*MRn1, *MRd1);
  
  s = "Num: ";
  s += MRn1->Integral();
  s1 = "Den: ";
  s1 += MRd1->Integral();
  std::cout << s << std::endl;
  std::cout << s1 << std::endl;
  
  MRn1->Divide(MRd1);
  MRn1->SetStats(0);
  MRn1->SetMarkerStyle(24);
  MRn1->SetMarkerSize(1.5);
  MRn1->SetMarkerColor(4);
  MRn1->SetTitle("");
  MRn1->SetXTitle("M_{R}");
  MRn1->GetYaxis()->SetRangeUser(0.0, 1.2);
  MRn1->Draw();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(22);
  tex->DrawLatex(0.8,0.8,"R^{2} > 0.3");
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.3,s);
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.2,s1);
  

  C1->SaveAs("SingleElectronPDTrigger_MR.png");
  C1->SaveAs("SingleElectronPDTrigger_MR.pdf");
  delete tex;

  pEff3->Draw();
  pEff3->SetMarkerStyle(24);
  pEff3->SetMarkerSize(1.5);
  pEff3->SetMarkerColor(4);
  pEff3->SetLineColor(4);
  pEff3->SetTitle(";M_{R};");
  pEff3->Draw();
  gPad->Update();
  pEff3->GetPaintedGraph()->GetXaxis()->SetTitle("M_{R}");
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(22);
  tex->DrawLatex(0.8,0.6,"R^{2} > 0.3");
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.3,s);
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.2,s1);
  delete tex;
  C1->SaveAs("EleEff_MR.pdf");
  C1->SaveAs("EleEff_MR.png");
  
  R2n1->Sumw2();
  R2d1->Sumw2();
  TEfficiency* pEff4 = new TEfficiency(*R2n1, *R2d1);
  s = "Num: ";
  s += R2n1->Integral();
  s1 = "Den: ";
  s1 += R2d1->Integral();
  std::cout << "Num: " << s << std::endl;
  std::cout << "Den: " << s1 << std::endl;
  R2n1->Divide(R2d1);
  R2n1->SetStats(0);
  R2n1->SetMarkerStyle(24);
  R2n1->SetXTitle("R^{2}");
  R2n1->SetMarkerSize(1.5);
  R2n1->SetMarkerColor(4);
  R2n1->SetTitle("");
  R2n1->GetYaxis()->SetRangeUser(0.0, 1.2);
  R2n1->Draw();
  

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(22);
  tex->DrawLatex(0.7,0.85,"M_{R} > 200 GeV");
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.3,s);
  delete tex;

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.2,s1);


  
  C1->SaveAs("SingleElectronPDTrigger_MR.png");
  C1->SaveAs("SingleElectronPDTrigger_MR.pdf");
  delete tex;

  pEff4->SetMarkerStyle(24);
  pEff4->SetMarkerSize(1.5);
  pEff4->SetMarkerColor(4);
  pEff4->SetLineColor(4);
  pEff4->SetTitle(";M_{R};");
  pEff4->Draw();
  gPad->Update();
  pEff4->GetPaintedGraph()->GetXaxis()->SetTitle("R^{2}");
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(22);
  tex->DrawLatex(0.8,0.6,"M_{R} > 200");
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.3,s);
  delete tex;
  
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(10);
  tex->DrawLatex(0.7,0.2,s1);
  delete tex;
  C1->SaveAs("EleEff_R2.pdf");
  C1->SaveAs("EleEff_R2.png");

  
  std::cout << "ll " << num1->GetBinContent(3,3) << std::endl;
  
  TRandom3* r = new TRandom3(0);
  
  for(int j = 0; j < 100; j++){
    double mr = r->Uniform(200.,2000.);
    double rsq = r->Uniform(0.2, 1.5);
    //std::cout << "MR:  " << mr << " R2: "  << rsq << "  find bin: " << FindBin(mr,rsq) << std::endl;
    //if( FindBin(mr,rsq) == 0 )std::cout << "bad" << std::endl;
  }
  
  
  TFile* f1 = new TFile("hlt_eff_SignleElePD.root", "RECREATE");
  //num->Write("h");
  num1->Write();
  R2n1->Write();
  MRn1->Write();
  pEff2->Write("Eff2d");
  pEff3->Write("EffMR");
  pEff4->Write("EffR2");
  f1->Close();
  
  return 0;
  
}


int scaleFactor(double MR, double RSQ){
  
  double MRhigh = 2000.;
  double MRlow_fixed = 200.;
  double MRlow = MRlow_fixed;
  int i = 0;
  while(1){//Binary Search!!
    i++;
    if( MR > (MRhigh - MRlow)/2. + MRlow){
      MRlow = (MRhigh - MRlow)/2. + MRlow;
      continue;
    }else if(  MR < (MRhigh - MRlow)/2. + MRlow ){
      MRhigh = (MRhigh - MRlow)/2. + MRlow;
      continue;
    }else{
      std::cout << "it: " <<  i << "  number is:  " << (MRhigh - MRlow)/2. + MRlow << std::endl;
      break;
    }
    
    
  }
  return 1;
};


double FindBin(double MR, double R2){
  
  int MRbin = -1;
  int R2bin = -1;
  
  const double R2A[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MRA[] = {200., 300., 400. ,500., 2000.};
  
  int Nbins = 4;
  
  for(int j = 0; j < Nbins; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
	R2bin = j+1;
	break;
      }    
    }
  }
  
  for(int j = 0; j < Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
	MRbin = j+1;
        break;
      }
    }
  }
  
  
  //return num1->GetBinContent( MRbin, R2bin );
  pEff2->GetEfficiency(pEff2->GetGlobalBin(MRbin, R2bin , 0));
};
