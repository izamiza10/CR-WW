// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
//#include "Rivet/Logging.hh"

#include <string.h>

/*

A RIVET routine for trying to understand what measurements would be useful for color reconnection
List of current included observables:

* Color reconnection studies

Written by Jennifer Roloff <jroloff2@gmail.com>

For more suggestions, please email Jennifer Roloff

*/

namespace Rivet {


  /// @brief MC_COLORRECONNECTION // placeholder number 1111111 for now
  class MC_COLORRECONNECTION : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_COLORRECONNECTION);

    // Given a list of bins, gives the bin which the value falls in
    // If it is below or above the last bin value, it will return -1
    int getBin(const double value, const std::vector<double> bins){
      for(unsigned int i=0; i<bins.size()-1; i++){
        if(value >= bins[i] && value < bins[i+1]) return i;
      }
      return -1;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Note: if we want to get quarks and gluons before hadronization, can use FinalPartons

      // The basic final-state projection: all final-state particles within the given eta acceptance
      // 10 degrees < theta < 170 degrees
      const FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 0.1*GeV);
      declare(fs, "particles");

      ChargedFinalState tracks(Cuts::pT > 0.2*GeV && Cuts::abseta < m_maxEtaTracks);
      declare(tracks, "tracks");

      // TODO: I'm not sure why the Durham algorithm seems to be failing, but we should figure that out...
      //FastJets jetDurham(fs, FastJets::Algo::DURHAM, m_jetRadius, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      FastJets jetDurham(fs, FastJets::Algo::ANTIKT, m_jetRadius, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetDurham, "DurhamJets");


      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Book histograms

      // One histogram for each of the two jet pairs
      Histo1DPtr tmpDijetMass1;
      Histo1DPtr tmpDijetMass2;
      std::string nameDijetMass1 = "dijet_mass_1";
      std::string nameDijetMass2 = "dijet_mass_2";
      book(_h_dijetMass1, nameDijetMass1,  40, 0, 120);
      book(_h_dijetMass2, nameDijetMass2,  40, 0, 120);

      // One histogram for tracking what particles are left after cuts
      std::string cutflow = "Cutflow";
      book(_h_cutflow, cutflow, 5, 0.5, 5.5);

      //ycut before event vetoes
      book(_h_ycutAll, "ycutAll", 60, 1e-6, 0.5);

      //Selected jets before cut
      book(_h_selectedJets, "selectedJets", 13, -0.5, 12.5); 

      //particle multiplicity before cuts
      book(_h_mult_before, "PreMultiplicity", 150, -0.5, 150.5);   

      //particle multiplicity after cuts
      book(_h_mult_after, "PostMultiplicity", 150, -0.5, 150.5);

      //Figure 5
      book(phi_rescaled, "ParticleFlow_PhiRescaled", 84, -0.2, 4.2);

      //Figure 6
      book(insideTotal, "Inside_W_Region", 24, -0.2, 1.2);
      book(OutsideTotal, "Outside_W_Region", 24, -0.2, 1.2);
      book(region_ratio, "Region_Ratio", 24, -0.2, 1.2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Most inspiration taken from https://arxiv.org/pdf/0704.0597
      // Some plots that would be good to add:
      // A cutflow of what event selections are leading to events being removed
      // A plot of the ycut of all durjet jets (before any cuts)
      // A plot of the number of selected jets (before any cuts)
      
      // A plot of the particle multiplicity (before any cuts)
      // A plot of the particle multiplicity (after all cuts)

      // A plot like Figure 6 (will probably be challenging), for all 4 regions

      // Things we don't need right now, but might want later
      //double collisionE = sqrtS();
      //const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();
      const Particles& particles = apply<FinalState>(event, "particles").particles();

      //All events
      //Use event.weight() maybe to normalize
      _h_cutflow->fill(1);

      // Fastjet analysis - select algorithm and parameters
      fastjet::Strategy               strategy = fastjet::Best;
      fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
      fastjet::JetDefinition*         jetDef = new fastjet::JetDefinition(fastjet::ee_kt_algorithm,
                                      recombScheme, strategy);

      // Run Fastjet algorithm
      fastjet::ClusterSequence clustSeq(particles, *jetDef);
      int nJetMin = 4;
      //clustSeq.exclusive_dmerge_max(nJetMin-1);

      const std::vector<fastjet::PseudoJet> durjet = clustSeq.exclusive_jets(nJetMin);

      _h_mult_before->fill(particles.size());

      std::vector<fastjet::PseudoJet> selectedJets;

      for(unsigned int i=0; i<durjet.size(); i++){
        PseudoJet jj, j1, j2;
        jj = durjet[i];
        if(jj.has_parents(j1, j2)){
          FourMomentum fv1 = FourVector(j1.E(), j1.px() , j1.py(), j1.pz());
          FourMomentum fv2 = FourVector(j2.E(), j2.px() , j2.py(), j2.pz());
          double ycut = 2 * j2.E()*j2.E() * (1- std::cos(fv1.angle(fv2))) / (jj.E()*jj.E());

          //Add ycut pre-veto
          _h_ycutAll->fill(ycut);

          if(ycut < 0.005) {
            selectedJets.push_back(jj);
          }
        }
      }

      //Histo of selected jets
      _h_selectedJets->fill(durjet.size());

      // Need exactly 4 jets
      if (selectedJets.size() != 4) vetoEvent;

      //Add event that survived 4-jet selection
      _h_cutflow->fill(2);

      // Two smallest angles < 100 degrees
      int minJetPair_1_jet1index = -1;
      int minJetPair_1_jet2index = -1;

      int minJetPair_2_jet1index = -1;
      int minJetPair_2_jet2index = -1;

      double minJetAngle1 = 1000;
      double minJetAngle2 = 1000;

      for(unsigned int i=0; i<selectedJets.size(); i++){
        for(unsigned int j=i; j<selectedJets.size(); j++){
          if(i==j) continue;
          FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
          FourMomentum fv2 = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          double angleInDegrees = fv1.angle(fv2) * 180 / 3.14;

          if(angleInDegrees < minJetAngle1){
            minJetPair_2_jet1index = minJetPair_1_jet1index;
            minJetPair_2_jet2index = minJetPair_1_jet2index;
            minJetAngle2 = minJetAngle1;

            minJetPair_1_jet1index = i;
            minJetPair_1_jet2index = j;
            minJetAngle1 = angleInDegrees;
          }
          else if(angleInDegrees < minJetAngle2){
            minJetPair_2_jet1index = i;
            minJetPair_2_jet2index = j;
            minJetAngle2 = angleInDegrees;
          }

        }
      }

      // Find two other angles that are between 100 and 140 degrees, and aren't adjacent
      //These become the jet pairs from the W
      std::vector<double> jetAngles;
      //std::vector<double> jetMasses;
      std::vector<int> jetIndices1;
      std::vector<int> jetIndices2;

      for(unsigned int i=0; i<selectedJets.size(); i++){
        for(unsigned int j=i+1; j<selectedJets.size(); j++){
          //The selected jet pairs aren't the smallest pairs
          if( (i == minJetPair_1_jet1index && j == minJetPair_1_jet2index) ||
              (j == minJetPair_1_jet1index && i == minJetPair_1_jet2index) ||
              (i == minJetPair_2_jet1index && j == minJetPair_2_jet2index) ||
              (j == minJetPair_2_jet1index && i == minJetPair_2_jet2index) ) continue;

          FourMomentum fv1 = FourVector(selectedJets[i].E(), selectedJets[i].px() , selectedJets[i].py(), selectedJets[i].pz());
          FourMomentum fv2 = FourVector(selectedJets[j].E(), selectedJets[j].px() , selectedJets[j].py(), selectedJets[j].pz());
          double angleInDegrees = fv1.angle(fv2) * 180 / 3.14;
          //fastjet::PseudoJet jetPair1 = selectedJets[i] + selectedJets[j];

          //Collect the 4 remaining jet combinations and their angles
          jetIndices1.push_back(i);
          jetIndices2.push_back(j);
          //jetMasses.push_back(jetPair1.m());
          jetAngles.push_back( angleInDegrees);
        }
      }

      //Apply cut that the smallest angles must be less than 100 degrees
      if(minJetAngle1 > 100) vetoEvent;
      if(minJetAngle2 > 100) vetoEvent;

      //Cut if the pairs of jets are too wide
      _h_cutflow->fill(3);

      // These two angles don't share a jet
      // Not sure if we should veto the event, or just find another pairing
      if( minJetPair_1_jet1index ==  minJetPair_2_jet1index || minJetPair_1_jet1index ==  minJetPair_2_jet2index ||
          minJetPair_1_jet2index ==  minJetPair_2_jet1index || minJetPair_1_jet2index ==  minJetPair_2_jet2index){
         vetoEvent;
      }

      //Cut if jets are not distinct pairs
      _h_cutflow->fill(4);

      bool isPass = false;
      int index1 = -1;
      int index2 = -1;
      int index3 = -1;
      int index4 = -1;

      //Angles between W jets should be between 100 and 140 degrees
      for(unsigned int i=0; i<jetAngles.size(); i++){
        if(jetAngles[i] < 100 || jetAngles[i] > 140) continue;
        for(unsigned int j=i+1; j<jetAngles.size(); j++){
          if(jetAngles[j] < 100 || jetAngles[j] > 140) continue;

          //Loop through all pairings until we find a non adjacent one
          if( jetIndices1[i] ==  jetIndices1[j] || jetIndices2[i] ==  jetIndices2[j] ||
              jetIndices1[i] ==  jetIndices2[j] || jetIndices2[i] ==  jetIndices1[j]) continue;

          isPass = true;
          //Indexes 1 and 2 should be the first pair, indexes 3 and 4 should be the second
          index1 = jetIndices1[i];
          index2 = jetIndices2[i];
          index3 = jetIndices1[j];
          index4 = jetIndices2[j];
        }
      }

      if(!isPass) vetoEvent;

      //cut if the jets don't fit expected geometry
      _h_cutflow->fill(5);

      //fastjet::PseudoJet correctJetPair1 = selectedJets[index1] + selectedJets[index2];
      //fastjet::PseudoJet correctJetPair2 = selectedJets[index3] + selectedJets[index4];

      //We're gonna try to ensure the jets are listed in clockwise order
      //The paper says that the interjet regions have the smallest angles, which we found.
      
      std::array<int,4> indices = {index1,index2,index3,index4};
      std::vector<int> ordered_indices;

      //Start by putting the first of the smallest angle region
      ordered_indices.push_back(minJetPair_1_jet1index);

      //loop until we find which index matches it, then add its pair
      for(unsigned int i=0; i<indices.size();i++){
        if(minJetPair_1_jet1index == indices[i]){
          if(i%2 == 0){
            ordered_indices.push_back(indices[i+1]);
            break;
          }
          else{
            ordered_indices.push_back(indices[i-1]);
            break;
          }
        }
      }

      //Add the next jet based on the smallest jet pairs
      if(ordered_indices[1] == minJetPair_2_jet1index)
        ordered_indices.push_back(minJetPair_2_jet2index);
      else
        ordered_indices.push_back(minJetPair_2_jet1index);
      
      //Add the last remaining jet, and the first jet again for looping
      ordered_indices.push_back(minJetPair_1_jet2index);
      ordered_indices.push_back(minJetPair_1_jet1index);
      
      //Ok now ordered_indices should contain all jets in order clockwise
      //Note: First two are from the same W, next 2 are the other W

      fastjet::PseudoJet correctJetPair1 = selectedJets[ordered_indices[0]] + selectedJets[ordered_indices[1]];
      fastjet::PseudoJet correctJetPair2 = selectedJets[ordered_indices[2]] + selectedJets[ordered_indices[3]];

      _h_dijetMass1->fill(correctJetPair1.m());
      _h_dijetMass2->fill(correctJetPair2.m());
      _h_mult_after->fill(particles.size());

      //A quest for Figure 6

      for (int i = 0; i < 4; ++i) {
        
        //First we will define a normal vector from a couple of 3 Vectors:
        Vector3 a(selectedJets[ordered_indices[i]].px(), selectedJets[ordered_indices[i]].py(), selectedJets[ordered_indices[i]].pz());
        Vector3 b(selectedJets[ordered_indices[i+1]].px(), selectedJets[ordered_indices[i+1]].py(), selectedJets[ordered_indices[i+1]].pz());

        Vector3 n = a.cross(b);
        double nmag = n.mod();
        Vector3 unit_n = n / nmag;

        //loop through all particles and see if its within the region
        for (const Particle& part : particles){
          Vector3 p(part.px(), part.py(), part.pz());
          Vector3 proj_p = p - unit_n * p.dot(unit_n);

          //Calculate angles with a as the reference jet
          double thetaRef = atan2(a.cross(b).dot(unit_n), a.dot(b));
          double theta_p = atan2(a.cross(proj_p).dot(unit_n), a.dot(proj_p));

          //if its in the region, calculate a phi rescaled and fill
          if ((0 < theta_p && theta_p < thetaRef) || (0 > theta_p && theta_p > thetaRef)){
            double thetaRescaled = theta_p/thetaRef;
            phi_rescaled->fill(thetaRescaled+i);

            //Check if inside or outside region and fill cooresponding histo
            if (i == 0 || i == 2){
              insideTotal->fill(thetaRescaled);
            }
            if (i == 1 || i == 3){
              OutsideTotal->fill(thetaRescaled);
            }
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      normalize(_h_dijetMass1);
      normalize(_h_dijetMass2);
      normalize(_h_ycutAll);
      normalize(_h_selectedJets);
      normalize(_h_mult_before);
      normalize(_h_mult_after);
      scale(phi_rescaled, 1.0 / sumOfWeights());
      scale(insideTotal, 1.0 / sumOfWeights());
      scale(OutsideTotal, 1.0 / sumOfWeights());

      //Find actual ratio for figure 6

      for (size_t i = 0; i < insideTotal->numBins(); ++i) {
          
          // Get the bin content
          double insideValue  = insideTotal->bin(i).sumW();
          double outsideValue = OutsideTotal->bin(i).sumW();

          // Compute the ratio safely
          double ratio = (outsideValue != 0) ? insideValue / outsideValue : 0;

          // Fill the ratio into the corresponding bin of region_ratio
          region_ratio->fillBin(i, ratio);
      } 
    }

    bool m_debug = false;

    // If you want to match LEP, you should use this value for the max track eta
    //const double m_maxEtaTracks = 1.738; // Using cos(theta) < 0.94, which this should correspond to

    // For these studies, currently studying the same eta cutoff for all paricles and for tracks
    // This choice is somewhat arbitrary, and we could probably go higher in eta with future detectors
    const double m_maxEtaTracks = 2.1;
    const double m_maxEta = 2.1;

    const double m_jetRadius = 0.7;

    //std::vector<Histo1DPtr> _h_dijetMasses;
    Histo1DPtr _h_dijetMass1;
    Histo1DPtr _h_dijetMass2;
    //Same thing, std::vector<Histo1DPtr> _h_cutflow;
    Histo1DPtr _h_cutflow;
    //ycut before event cuts
    Histo1DPtr _h_ycutAll;
    //# of selected jets before cut
    Histo1DPtr _h_selectedJets;
    //particle multiplicity before cuts
    Histo1DPtr _h_mult_before;
    //particle multiplicity after cuts
    Histo1DPtr _h_mult_after;

    //Figure 6, particle flow per event in each region
    Histo1DPtr phi_rescaled;

    //Figure 5, ratio between particle count inside/outside
    Histo1DPtr insideTotal;
    Histo1DPtr OutsideTotal;
    Histo1DPtr region_ratio;
  };


  RIVET_DECLARE_PLUGIN(MC_COLORRECONNECTION);

}


