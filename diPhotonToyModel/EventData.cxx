#include "EventData.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

// generate the events

// parameters

int NTRACKS = 50; 
// parameter for pt distribution
double PtAvg = 0.6;   
int NTYPES = 3;
double masses[] = { 0, 0.0005, 0.1350}; //photon, election+/-, pi0
int charge [] = { 0, 1, 0}; 
double fractions[] = { 0.1, 0.1, 0.8};  
double sigmax = 10.E-6; 
double sigmay = 10.E-6; 
double sigmaz = 5.; 

using namespace ROOT::Math;

ROOT::Math::PtEtaPhiMVector SmearVector(const ROOT::Math::XYZTVector & v) { 
   double x = v.X()*(1. + gRandom->Gaus(0, 0.05) );
   double y = v.Y()*(1. + gRandom->Gaus(0, 0.05) );
   double z = v.Z()*(1. + gRandom->Gaus(0, 0.05) );
   ROOT::Math::PxPyPzMVector tmp(x,y,z,v.M() );
   return PtEtaPhiMVector(tmp);
}


void EventData::Generate()  { 


   // get expected value for each type
   for (int i = 0; i < NTYPES; ++i) { 
      double nexp = fractions[i] * NTRACKS; 

      int np = gRandom->Poisson(nexp); 
      for (int j = 0; j < np; ++j) { 
         Particle p; 
         p.fPosition = XYZVector( gRandom->Gaus(0,sigmax), gRandom->Gaus(0, sigmay), gRandom->Gaus(0, sigmaz) ); 
         double pt = gRandom->Exp(PtAvg); 
         double eta = gRandom->Uniform(-3,3);
         double phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi() ); 
         double mass = masses[i];
         p.fVector = PtEtaPhiMVector(pt, eta, phi, mass); 
         p.fType = i; 
         p.fCharge = charge[i];
         if (p.fCharge) { 
            int tmp = gRandom->Integer(2); 
            if (tmp == 0) p.fCharge = -1; 
         }
         // special case for decays
         if (i == 2 ) { 
            // pi0 to two photons
            TGenPhaseSpace evt; 
            double m[2] = {0,0};
            TLorentzVector W( p.fVector.X(), p.fVector.Y(), p.fVector.Z(), p.fVector.E() );
            evt.SetDecay(W, 2, m);
            evt.Generate();
            TLorentzVector * v1 = evt.GetDecay(0); 
            TLorentzVector * v2 = evt.GetDecay(1); 
            Particle p1; 
            Particle p2; 
            p1.fPosition = p.fPosition; 
            p2.fPosition = p.fPosition; 
            p1.fCharge = 0;
            p2.fCharge = 0;
            p1.fType = 0; 
            p2.fType = 0;

            p1.fVector = SmearVector(XYZTVector(v1->X(), v1->Y(), v1->Z(), v1->E() )); 
            p2.fVector = SmearVector(XYZTVector(v2->X(), v2->Y(), v2->Z(), v2->E() )); 
            AddParticle(p1);
            AddParticle(p2);
         }
         else
            AddParticle(p);

      }
   }
}
