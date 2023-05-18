//Running the code:
//[tuos@gw341 test]$ root -l
//root [0] .L CreateEventTree.C+
//root [1] .x CreateEventTree.C
//
// Prefer compiled:
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>

#include "EventData.h"
#include "EventData.cxx"

void CreateEventTree(Long64_t numEvents = 100000, int splitlevel = 99) {
   TFile* f = new TFile("./eventdata_100k.root", "RECREATE");
   TTree* tree = new TTree("tree", "Tutorial tree");

   EventData* event = new EventData();
   tree->Branch("event", &event,32000,splitlevel);

   Particle p;

   for (Long64_t i = 0; i < numEvents; ++i) {
      if(i%1000==0) cout<<"i = "<<i<<" ,  numEvents="<<numEvents<<endl;
      event->Clear();
      event->Generate();
      event->SetSize();

      tree->Fill();

   }

   tree->Print();
   tree->Write();
   delete f;
}
