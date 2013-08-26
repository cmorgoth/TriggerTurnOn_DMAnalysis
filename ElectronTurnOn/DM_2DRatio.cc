#include "DM_2DRatio.hh"

TH2F* RatioPlotsTwoD(TH2F* h1, TH2F* h2, TString fname = "fname", TString hname = "histo_name"){
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  TH2F* Ratio2D = new TH2F( hname, hname, BaseDM::MR_Bins, BaseDM::MR_BinArr, BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
  Ratio2D->Divide( h1, h2, 1, 1, "");
  
  Ratio2D->SetStats(0);
  Ratio2D->SetTitle("");
  Ratio2D->GetXaxis()->SetRangeUser(200.0, 3500.);
  Ratio2D->GetYaxis()->SetRangeUser(0.3, 2.50);
  Ratio2D->SetXTitle("M_{R}");
  Ratio2D->SetYTitle("R^{2}");
  Ratio2D->Draw("colz");
  
  float binV, binErr, Xpos, Ypos;
  TLatex* Lval;
  TLatex* Lerr;
  char val[50];
  char err[50];
  binV = 666.;
  /////////////////////////////////////////////////////
  ////////////////////2D PLOTS with numbers///////////
  ////////////////////////////////////////////////////
  for(int i = 1; i <= Ratio2D->GetNbinsX(); i++){
    for(int j = 1; j <= Ratio2D->GetNbinsY(); j++){
      binV =  Ratio2D->GetBinContent(i,j);
      binErr = Ratio2D->GetBinError(i,j);
      
      Xpos =  Ratio2D->GetXaxis()->GetBinCenter(i);
      Ypos =  Ratio2D->GetYaxis()->GetBinCenter(j);
      
      sprintf(val,"%.3f", binV);
      sprintf(err," \\pm %.3f", binErr);
      //std::cout << "val: " << val << err << "Ypos: " << Ypos << std::endl;
      if( (binV != 0) || ( binErr != 0) ){
	Lval = new TLatex(Xpos-70,Ypos+0.01, "HELLO!!!" );
	Lerr = new TLatex(Xpos-70,Ypos-0.03, err );
	Lval->SetTextSize(.015);
	Lval->Draw("same");
	Lval->AppendPad();
	Lval->DrawLatex(Xpos -12. ,Ypos+0.01, val );
	Lval->DrawLatex(Xpos -12. ,Ypos-0.013, err );
	Lerr->Draw("same");
	
	delete Lval;
	delete Lerr;
	//std::cout << printf( "%.4f", binV ) << std::endl;
      }
     
    }
    
  }
  
  C1->SaveAs( fname + ".pdf" );
  C1->SaveAs( fname + ".png" );
  
  delete C1;
  
  return Ratio2D;
  

};
