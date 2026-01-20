#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include <algorithm>
#include <numeric>

namespace Rivet {

class MY_ANALYSIS_WW_truthDecay : public Analysis {
public:
  RIVET_DEFAULT_ANALYSIS_CTOR(MY_ANALYSIS_WW_truthDecay);

  void init() override {

    book(_h_deltaThetaWp, "quark_deltaTheta_Wp", 50, 0.0, M_PI);
    book(_h_deltaPhiWp, "quark_deltaPhi_Wp", 50, -M_PI, M_PI);
    book(_h_deltaRWp, "quark_deltaR_Wp", 50, 0.0, 6.0);

    book(_h_deltaThetaWm, "quark_deltaTheta_Wm", 50, 0.0, M_PI);
    book(_h_deltaPhiWm, "quark_deltaPhi_Wm", 50, -M_PI, M_PI);
    book(_h_deltaRWm, "quark_deltaR_Wm", 50, 0.0, 6.0);

    book(_h_deltaThetaInt1, "interjet_deltaTheta_q1q3", 50, 0.0, M_PI);
    book(_h_deltaPhiInt1, "interjet_deltaPhi_q1q3", 50, -M_PI, M_PI);
    book(_h_deltaRInt1, "interjet_deltaR_q1q3", 50, 0.0, 6.0);

    book(_h_deltaThetaInt2, "interjet_deltaTheta_q2q4", 50, 0.0, M_PI);
    book(_h_deltaPhiInt2, "interjet_deltaPhi_q2q4", 50, -M_PI, M_PI);
    book(_h_deltaRInt2, "interjet_deltaR_q2q4", 50, 0.0, 6.0);

    book(_h_angle_Wp, "openingAngle_Wp", 50, 0.0, 180.0);
    book(_h_angle_Wm, "openingAngle_Wm", 50, 0.0, 180.0);
    book(_h_angle_Int1, "openingAngle_q1q3", 50, 0.0, 180.0);
    book(_h_angle_Int2, "openingAngle_q2q4", 50, 0.0, 180.0);

    book(_h_massWp, "Wplus_mass", 50, 0.0, 200.0);
    book(_h_massWm, "Wminus_mass", 50, 0.0, 200.0);
  }

  void analyze(const Event &event) override {

    const HepMC3::GenEvent &ge = *event.genEvent();
    //  Find hard-scatter W bosons
    HepMC3::ConstGenParticlePtr Wp = nullptr;
    HepMC3::ConstGenParticlePtr Wm = nullptr;

    for (const auto &p : ge.particles()) {
      if (abs(p->pid()) != 24)
        continue;
      if (!p->production_vertex())
        continue;
      bool prod_ee = true;
      for (const auto &parent : p->production_vertex()->particles_in()) {
        if (abs(parent->pid()) != 11) {
          prod_ee = false;
          break;
        }
      }
      if (!prod_ee)
        continue;

      if (p->pid() == 24)
        Wp = p;
      if (p->pid() == -24)
        Wm = p;
    }
    if (!Wp || !Wm)
      return;

    //  Get quarks from W decays
    std::vector<Particle> WpQuarks, WmQuarks;

    if (Wp->end_vertex()) {
      for (const auto &q : Wp->end_vertex()->particles_out())
        if (PID::isQuark(q->pid()))
          WpQuarks.push_back(Particle(q));
    }
    if (Wm->end_vertex()) {
      for (const auto &q : Wm->end_vertex()->particles_out())
        if (PID::isQuark(q->pid()))
          WmQuarks.push_back(Particle(q));
    }
    if (WpQuarks.size() != 2 || WmQuarks.size() != 2)
      return;

    //  Order Wp quarks by pT
    if (WpQuarks[1].pT() > WpQuarks[0].pT())
      std::swap(WpQuarks[0], WpQuarks[1]);

    Particle q1 = WpQuarks[0];
    Particle q2 = WpQuarks[1];

    // Pick q3 based on mind dR to q1
    double dR0 = deltaR(q1.momentum(), WmQuarks[0].momentum());
    double dR1 = deltaR(q1.momentum(), WmQuarks[1].momentum());

    Particle q3 = (dR0 < dR1 ? WmQuarks[0] : WmQuarks[1]);
    Particle q4 = (dR0 < dR1 ? WmQuarks[1] : WmQuarks[0]);

    // Setup 4-momentum
    FourMomentum fv_q1(q1.E(), q1.px(), q1.py(), q1.pz());
    FourMomentum fv_q2(q2.E(), q2.px(), q2.py(), q2.pz());
    FourMomentum fv_q3(q3.E(), q3.px(), q3.py(), q3.pz());
    FourMomentum fv_q4(q4.E(), q4.px(), q4.py(), q4.pz());

    double angleWpDeg = fv_q1.angle(fv_q2) * 180.0 / M_PI;
    double angleWmDeg = fv_q3.angle(fv_q4) * 180.0 / M_PI;

    // Only fill if opening angle is 100-140 degrees, missing the cut on
    // interjet region
    if (angleWpDeg > 100.0 && angleWpDeg < 140.0 && angleWmDeg > 100.0 &&
        angleWmDeg < 140.0) {
      _h_angle_Wp->fill(angleWpDeg);
      _h_massWp->fill((fv_q1 + fv_q2).mass());
      fillQuarkHistograms(q1, q2, _h_deltaThetaWp, _h_deltaPhiWp, _h_deltaRWp);

      _h_angle_Wm->fill(angleWmDeg);
      _h_massWm->fill((fv_q3 + fv_q4).mass());
      fillQuarkHistograms(q3, q4, _h_deltaThetaWm, _h_deltaPhiWm, _h_deltaRWm);

      // Interjet angles
      _h_angle_Int1->fill(fv_q1.angle(fv_q3) * 180.0 / M_PI);
      _h_angle_Int2->fill(fv_q2.angle(fv_q4) * 180.0 / M_PI);
      fillQuarkHistograms(q1, q3, _h_deltaThetaInt1, _h_deltaPhiInt1,
                          _h_deltaRInt1);
      fillQuarkHistograms(q2, q4, _h_deltaThetaInt2, _h_deltaPhiInt2,
                          _h_deltaRInt2);
    }
  }

  void finalize() override {
    normalize(_h_deltaThetaWp);
    normalize(_h_deltaPhiWp);
    normalize(_h_deltaRWp);
    normalize(_h_deltaThetaWm);
    normalize(_h_deltaPhiWm);
    normalize(_h_deltaRWm);
    normalize(_h_deltaThetaInt1);
    normalize(_h_deltaPhiInt1);
    normalize(_h_deltaRInt1);
    normalize(_h_deltaThetaInt2);
    normalize(_h_deltaPhiInt2);
    normalize(_h_deltaRInt2);
    normalize(_h_angle_Wp);
    normalize(_h_angle_Wm);
    normalize(_h_angle_Int1);
    normalize(_h_angle_Int2);
    normalize(_h_massWp);
    normalize(_h_massWm);
  }

private:
  void fillQuarkHistograms(const Particle &a, const Particle &b,
                           Histo1DPtr &hTheta, Histo1DPtr &hPhi,
                           Histo1DPtr &hR) {

    double deltaPhi = a.phi() - b.phi();
    while (deltaPhi > M_PI)
      deltaPhi -= 2 * M_PI;
    while (deltaPhi < -M_PI)
      deltaPhi += 2 * M_PI;

    double deltaEta = a.eta() - b.eta();
    double deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

    hPhi->fill(deltaPhi);
    hTheta->fill(std::abs(a.theta() - b.theta()));
    hR->fill(deltaR);
  }

  Histo1DPtr _h_deltaThetaWp, _h_deltaPhiWp, _h_deltaRWp;
  Histo1DPtr _h_deltaThetaWm, _h_deltaPhiWm, _h_deltaRWm;
  Histo1DPtr _h_deltaThetaInt1, _h_deltaPhiInt1, _h_deltaRInt1;
  Histo1DPtr _h_deltaThetaInt2, _h_deltaPhiInt2, _h_deltaRInt2;
  Histo1DPtr _h_angle_Wp, _h_angle_Wm, _h_angle_Int1, _h_angle_Int2;
  Histo1DPtr _h_massWp, _h_massWm;
};

RIVET_DECLARE_PLUGIN(MY_ANALYSIS_WW_truthDecay);

} // namespace Rivet
