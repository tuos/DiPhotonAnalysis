// program to test treereader functionality
//Running the code
//[tuos@gw341 test]$ root -l
//root [0] .L EventDataReader.C+
//root [1] .x EventDataReader.C 
//

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/LorentzVector.h"
//#include "Math/Vector4d.h"
#include "Math/Vector4D.h"
#include <iostream>

#include "TStopwatch.h"

#include "TCanvas.h"
#include "EventData.h"

void EventDataReader() {

   TH1::SetDefaultSumw2();
   TH1 * h1 = new TH1F("h1","inv mass pi0",90,0,0.3);
   TH1 * h2 = new TH1F("h2","inv mass pi0 with photon conversions",90,0,0.3);
   TH1 * h3 = new TH1F("h3","inv mass photon conversions only",90,0,0.3);

   TFile *myFile = TFile::Open("./eventdata_photonConversion_100k.root");

   TTreeReader myReader("tree", myFile);

   TTreeReaderValue<std::vector< Particle> > particles_value(myReader, "fParticles");


   std::cout << myReader.GetEntries(1) << std::endl;
   Long64_t totalEvent = myReader.GetEntries(1);

   TStopwatch w; 
   w. Start(); 
   int iEvents = 0;
   //while (myReader.Next()) {
   while (myReader.Next() && iEvents<10000) {
      std::vector< Particle>  & fParticles = *particles_value; 
      int npart = fParticles.size(); 

      iEvents++;
      if(iEvents%1000==0) cout<<"iEvents = "<<iEvents<<" ,  totalEvents="<<totalEvent<<endl;

      for (int i = 0; i< npart; ++i) { 
         for (int j = i+1; j< npart; ++j) { 
            if (fParticles[j].fType == 0 && fParticles[i].fType == 0) 
               h1->Fill( (fParticles[i].fVector + fParticles[j].fVector).M() );
            if ( (fParticles[j].fType == 0 && fParticles[i].fType == 0) || 
                 ( (fParticles[i].fType==0 && fParticles[j].fType  == 1) &&
                   (fParticles[i].fVector.Dot(fParticles[j].fVector) < 0.001) 
                 ) ) 
               h2->Fill( (fParticles[i].fVector + fParticles[j].fVector).M() );
            if ( (fParticles[i].fType==0 && fParticles[j].fType  == 1) && 
                 (fParticles[i].fVector.Dot(fParticles[j].fVector) < 0.001) ) 
               h3->Fill( (fParticles[i].fVector + fParticles[j].fVector).M() );

         }
      }

      
   }

   std::cout << "Time to read the tree:  ";
   w.Print();
   TCanvas * c1 = new TCanvas(); 
   c1->Divide(2,1);
   c1->cd(1);
   h1->SetLineColor(1);
   h1->SetMarkerColor(1);
   h1->SetMarkerStyle(20);
   h1->SetMarkerSize(1);
   h1->Draw("pe");
   h2->SetLineColor(2);
   h2->SetMarkerColor(2);
   h2->SetMarkerStyle(24);
   h2->SetMarkerSize(1);
   h2->Draw("pesame");
   c1->cd(2);
   h3->SetLineColor(4);
   h3->SetMarkerColor(4);
   h3->SetMarkerStyle(25);
   h3->SetMarkerSize(1);
   h3->Draw("pe");

   c1->SaveAs("fig_pi0_diphoton_conversions.png");

}


