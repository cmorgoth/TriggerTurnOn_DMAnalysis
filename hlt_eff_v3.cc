#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include <vector>
#include <iostream>

int scaleFactor(double MR, double RSQ);
double FindBin(double MR, double R2);

//TFile* f;
TH2F* num1;
TEfficiency* pEff2;

int main(){

  gROOT->Reset();
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  const double RSQ_BinArr[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MR_BinArr[] = {200., 300., 400. ,500., 2000.};
    
  TFile* f = new TFile("MuHad_ABCD_PT80Total.root");
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
  
  num1 = new TH2F("num1","num1", 4, MR_BinArr, 4, RSQ_BinArr);
  TH2F* den1 = new TH2F("den1","den1", 4, MR_BinArr, 4, RSQ_BinArr);
  
  TH1F* MRn1 = new TH1F("MRn1","Muhad-Turn-On-MR", 4, MR_BinArr);
  TH1F* MRd1 = new TH1F("MRd1","Muhad-Turn-On-MR-Den", 4, MR_BinArr);

  TH1F* R2n1 = new TH1F("R2n1","Muhad-Turn-On-R2", 4, RSQ_BinArr);
  TH1F* R2d1 = new TH1F("R2d1","Muhad-Turn-On-R2-Den", 4, RSQ_BinArr);

  
  
  std::cout << "Nentries: " << t->GetEntries() << std::endl;
  
  for(int i = 0; i < t->GetEntries(); i++ ){
    t->GetEntry(i);
    //if( RSQ[2] >= 0.35 && MR[2] >= 150. && HLT_Razor_prescaled == 1 && BOX >= 1){
    if( RSQ[2] >= 0.35 && MR[2] >= 200. && BOX >= 1 ){
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
  
  //t->Draw("RSQ[2]:MR[2]>>tmp1(50,200,2000,50, 0.5, 1.5)", "RSQ[2]>=0.5 && MR[2]>=200. && passedHLT==1", "goff");
  //num = (TH2F*)gDirectory->Get("tmp1");
  //t->Draw("RSQ[2]:MR[2]>>tmp2(50,200,2000,50, 0.5, 1.5)", "RSQ[2]>=0.5 && MR[2]>=200.", "goff");
  //TH2F* den = (TH2F*)gDirectory->Get("tmp2");
 
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
  C1->SaveAs("mu_trigger_2d_mr200_ABCD_PT80v2Muon.png");
  C1->SaveAs("mu_trigger_2d_mr200_ABCD_PT80v2Muon.pdf");

  pEff2->Draw("colz");
  C1->SaveAs("Eff_Tets_2d.pdf");
  C1->SaveAs("Eff_Tets_2d.png");

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
  std::cout << "Den BinContent(4): " << MRn1->GetBinContent(4) << std::endl;
  std::cout << "Den BinContent(4): " << MRd1->GetBinContent(4) << std::endl;

  MRn1->Divide(MRd1);
  MRn1->SetStats(0);
  MRn1->SetMarkerStyle(24);
  MRn1->SetMarkerSize(1.5);
  MRn1->SetMarkerColor(4);
  MRn1->SetTitle("");
  MRn1->SetXTitle("M_{R}");
  MRn1->GetYaxis()->SetRangeUser(0.5, 1.2);
  MRn1->Draw();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSizePixels(22);
  tex->DrawLatex(0.8,0.8,"R^{2} > 0.35");
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
  

  C1->SaveAs("mu_trigger_MR_mr200_ABCD_PT80v2Muon.png");
  C1->SaveAs("mu_trigger_MR_mr200_ABCD_PT80v2Muon.pdf");
  delete tex;

  pEff3->Draw();
  C1->SaveAs("Eff_Tets_MR.pdf");
  C1->SaveAs("Eff_Tets_MR.png");
  
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
  R2n1->GetYaxis()->SetRangeUser(0.5, 1.2);
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


  
  C1->SaveAs("mu_trigger_R2_mr200_ABCD_PT80v2Muon.png");
  C1->SaveAs("mu_trigger_R2_mr200_ABCD_PT80v2Muon.pdf");
  delete tex;

  pEff4->Draw();
  C1->SaveAs("Eff_Tets_R2.pdf");
  C1->SaveAs("Eff_Tets_R2.png");

  
  std::cout << "ll " << num1->GetBinContent(3,3) << std::endl;
  
  TRandom3* r = new TRandom3(0);
  
  for(int j = 0; j < 100; j++){
    double mr = r->Uniform(200.,2000.);
    double rsq = r->Uniform(0.5, 1.5);
    std::cout << "MR:  " << mr << " R2: "  << rsq << "  find bin: " << FindBin(mr,rsq) << std::endl;
    if( FindBin(mr,rsq) == 0 )std::cout << "bad" << std::endl;
  }
  
  
  TFile* f1 = new TFile("hlt_eff_mr200_MoreBin_ABCD_PT80v2Muon.root", "RECREATE");
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
