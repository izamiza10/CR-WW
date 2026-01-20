// main302.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: colour reconnection; e+e- events;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Plot the colour reconnection rate in W+W- hadronic events
// as a function of collision energy.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "Pythia8Plugins/InputParser.h"

using namespace Pythia8;

//==========================================================================

int main(int argc, char** argv){

  // Set up command line options.
  InputParser ip("Illustrates how to do matching and merging.",
    {"./main302 -c main302ckkwl.cmnd",
        "./main302 -c main302powheg.cmnd",});
  ip.require("c", "Use this user-written command file.", {"-cmnd"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Input file.
  string cmndFile = ip.get<string>("c");


  // Number of events.
  int nEvent = 1000000;

  // Histograms.
  Hist recrateSK1("reconnection rate SK1", 40, 0., 400.);
  Hist recrateSK2("reconnection rate SK2", 40, 0., 400.);

  //Adding an array of energy values
  vector<int> energies = {183, 189, 200, 206, 207};

  // Loop over reconnection model and CM energy.
  //for (int iRec = 0; iRec < 5; ++iRec) {
    //Hist& recrate = (iRec == 1) ? recrateSK1 : recrateSK2;
    //for (int eCM : energies) {

      // Quiet generator creation. Shorthand for the event record.
      Pythia pythia ("../share/Pythia8/xmldoc", false);
      pythia.readFile(cmndFile);
      Event& event = pythia.event;

      // Incoming beams and energy. Optionally no ISR.
      //pythia.readString("Beams:idA = -11");
      //pythia.readString("Beams:idB = 11");
      //pythia.settings.parm("Beams:eCM", eCM);
      //pythia.readString("PDF:lepton = off");

      // Hard process, with hadronic W decays.
      //pythia.readString("WeakDoubleBoson:ffbar2WW = on");
      //pythia.readString("24:onMode = off");
      //pythia.readString("24:onIfAny = 1 2 3 4 5");

      // Switch on CR.
      //pythia.readString("ColourReconnection:reconnect = on");
      //pythia.settings.mode("ColourReconnection:mode", iRec);
      //pythia.readString("ColourReconnection:forceResonance = on");

      // Reduce printout. Switch off hadronization and decays.
      //pythia.readString("Print:quiet = on");
      pythia.readString("HadronLevel:Hadronize = on");
      pythia.readString("HadronLevel:Decay = on");

      // If Pythia fails to initialize, exit with error.
      if (!pythia.init()) return 1;

      // [ADDED FOR HEPMC3 OUTPUT]
      // Create converter + output file (name encodes mode and energy)
      //std::string filename = "WW_CR_mode" + std::to_string(iRec)
      //                     + "_E" + std::to_string(eCM) + ".hepmc";
      std::string filename = cmndFile.replace(cmndFile.find(".cmnd"), 5, "") + ".hepmc";
      HepMC3::Pythia8ToHepMC3 toHepMC;
      HepMC3::WriterAscii writer(filename);

      int nRecon = 0;

      // Begin event loop.
      for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;

        // [ADDED FOR HEPMC3 OUTPUT]
        // Convert event and write to HepMC3 file
        HepMC3::GenEvent hepmcEvent;
        toHepMC.fill_next_event(pythia, &hepmcEvent);
        writer.write_event(hepmcEvent);

      }

      // [ADDED FOR HEPMC3 OUTPUT]
      writer.close();
      
      // End of loops. Fill reconnection rate.
      //recrate.fill( eCM, double(nRecon) / double(nEvent));
  //  }
  //}

  // Done.
  return 0;
}
