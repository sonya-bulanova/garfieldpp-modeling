#include <iostream>
#include <fstream>
#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

#include "TFile.h"
#include "TH2D.h"

//#include  "Garfield/Microscopic.hh""
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  MediumMagboltz gas;
  gas.LoadGasFile("Ar_70_co2_30_1atm_293.15K_B0.gas");
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << "GARFIELD_INSTALL variable not set.\n";
    return 1;
  }
  gas.LoadIonMobility("/home/bulanova.sa/University/garfieldpp/Data/IonMobility_Ar+_Ar.txt");

 const double rPenning = 0;
  // Mean distance from the point of excitation.
  const double lambdaPenning = 0.;
  //gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Make a component with analytic electric field.

  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Wire radius [cm]
  const double rWire = 15.e-4;
  // Outer radius of the xtube [cm]
  const double rTube = 0.50;
  // Voltages
  const double vWire = 1750.;
  const double vTube = 0.0;
  // Add the wire in the centre.
  cmp.AddWire(0, 0, 2 * rWire, vWire, "s");
  // Add the tube.
  cmp.AddTube(rTube, vTube, 0, "t");
  // Request calculation of the weighting field.
  cmp.AddReadout("s");
cmp.SetMagneticField(0,0,0);
  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "s");
  // Set the signal time window.
  const unsigned int nbins = 2000;
  const double tmin = 0;
  const double tmax = 250.;
  const double tstep = (tmax - tmin) / nbins;
 sensor.SetTimeWindow(tmin, tstep, nbins);

  TH1D hsignal("Name","Tittle", nbins, tmin, tmax);


  // Set up Heed.
  TrackHeed track;
  track.SetParticle("muon");
  track.SetEnergy(1.e9);
  track.SetSensor(&sensor);

  DriftLineRKF drift;

drift.SetSensor(&sensor);
drift.SetGainFluctuationsPolya(1,45000, false);
drift.EnableIonTail();
drift.EnableSignalCalculation();
  TCanvas* cD = nullptr;
  ViewCell cellView;
  ViewDrift driftView;
  constexpr bool plotDrift = false;
  if (plotDrift) {
    cD = new TCanvas("cD", "", 600, 600);
    cellView.SetCanvas(cD);
    cellView.SetComponent(&cmp);
    driftView.SetCanvas(cD);
    track.EnablePlotting(&driftView);
      drift.EnablePlotting(&driftView);
  }

  TCanvas* cS = nullptr;
  ViewSignal signalView;
  constexpr bool plotSignal = false;
  if (plotSignal) {
    cS = new TCanvas("cS", "", 600, 600);
    signalView.SetCanvas(cS);
    signalView.SetSensor(&sensor);
    signalView.SetLabelY("signal [µa]");

  }

   std::ofstream outfile;

  const double rTrack = 0.01; //cm
  const double x0 = rTrack;
  const double y0 = -sqrt(rTube * rTube - rTrack * rTrack)+0.001;
  const unsigned int nTracks = 5000;
  for (unsigned int j = 0; j < nTracks; ++j) {
    sensor.ClearSignal();
   // std::vector<std::array<double, 4> > electrons;
  track.NewTrack(x0, y0, 0, 0, 0, 1, 0);
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
    int nc = 0;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
     for (int k = 0; k < nc; ++k) {
        double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
        double dx = 0., dy = 0., dz = 0.;
        track.GetElectron(k, xe, ye, ze, te, ee, dx, dy, dz);
        drift.DriftElectron(xe, ye, ze, te);

        double ne = 0., ni = 0.;
   drift.GetAvalancheSize(ne, ni);
   //std::cout<<"NE = "<<ne<<std::endl;
    }
  }
  std::string outname="signal_out/Ar_70_co2_30_1750V_B0_1atm_0_1mm_bin2000_" + std::to_string(j)+".sig";

  outfile.open(outname, std::ios::out);

   for (unsigned int i = 0; i < nbins; ++i) {
   const double t = (i + tmin) * tstep;
   const double f = sensor.GetSignal("s", i);
   outfile << t/1000000000.0<< " " << f/1000000.0<< std::endl;
   }
  outfile.close();
}

 TCanvas s(" "," ",600,600);
 hsignal.Draw();
 std::cout<<"done"<<std::endl;
  app.Run(kTRUE);

}
