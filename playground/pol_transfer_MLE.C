/////////////////////////////////////////////////////////////////////
///////////////////////  pol_transfer_ML.C   ////////////////////////
/////////////////////////////////////////////////////////////////////
// Maximum-likelihood extraction of target polarization components //
/////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMinuit.h"
#include "TMatrixD.h"

#include <vector>
#include <cmath>
#include <iostream>

// ----------------------------------------------------------------------------
// 1) Event coefficients stored for the MLE
// ----------------------------------------------------------------------------
struct EventCoeff {
  double lambda0;    // false-asymmetry Fourier series evaluated at phi_i
  double lam[3];     // lambda_j coefficients multiplying (Pt,Pn,Pl)
  double weight;     // event weight (usually 1.0)
};

std::vector<EventCoeff> g_events;

// Global beam polarization 
double g_beam_polarization = 0.85;  // TODO: set for per run if possible

// ----------------------------------------------------------------------------
// 2) Physics helper stubs 
// ----------------------------------------------------------------------------

// Analyzing power A_y(p,theta)
double GetAnalyzingPower(double p, double th)
{
  return 1.0; // for now
}

// Compute spin-transport elements 
void ComputeSpinTransportForEvent(
    double Sx[3],  // out: S_xt, S_xn, S_xl
    double Sy[3]   // out: S_yt, S_yn, S_yl
)
{
  // TODO: replace with COSY-based spin transport using event kinematics
  Sx[0] = 0.0; Sx[1] = 0.0; Sx[2] = 0.0;
  Sy[0] = 1.0; Sy[1] = 0.0; Sy[2] = 0.0;
}

// False-asymmetry Fourier coefficients
void GetFalseAsymmetryCoeffs(double phi,
                            double& a1, double& b1,
                            double& a2, double& b2)
{
  // for now using zero false asymmetry.
  a1 = 0.0; b1 = 0.0;
  a2 = 0.0; b2 = 0.0;
}

// ----------------------------------------------------------------------------
// 3) Build the event list g_events from a TTree
// ----------------------------------------------------------------------------
void BuildEventList(TTree* T)
{
  g_events.clear();

  // Use TTreeReader for safer branch access
  TTreeReader Tl(T);

  TTreeReaderValue<double> helicity_branch(Tl, "heli");  
  TTreeReaderValue<double> phi_branch     (Tl, "phsc");  

  Long64_t nentries = T->GetEntries();
  std::cout << "[BuildEventList] Processing " << nentries << " events\n";

  while (Tl.Next()) {

    // Helicity eps = ±1
    // const double hraw = *helicity_branch;
    // const int eps = (hraw >= 0.0) ? +1 : -1;
    const int eps = *helicity_branch;
    const double h = g_beam_polarization;
    const double phi = *phi_branch;

    // Get Analyzing Power with  FPP kinematics 
    const double p_fpp  = 0.0;
    const double th_fpp = 0.0;

    const double Ay = GetAnalyzingPower(p_fpp, th_fpp);

    // Spin transport
    double Sx[3], Sy[3];
    ComputeSpinTransportForEvent(Sx, Sy);

    // Build raw lambda_j:
    //   lambda_j(raw) = Ay * ( Syj*cos(phi) - Sxj*sin(phi) )
    const double c = std::cos(phi);
    const double s = std::sin(phi);

    double lambda_raw[3];
    for (int j = 0; j < 3; j++) {
      lambda_raw[j] = Ay * ( Sy[j]*c - Sx[j]*s );
    }

    // Apply helicity dependence:
    const double hs = h * double(eps);

    double lambda_eff[3];
    lambda_eff[0] = hs * lambda_raw[0];   // Pt
    lambda_eff[1] =      lambda_raw[1];   // Pn (NO h*eps)
    lambda_eff[2] = hs * lambda_raw[2];   // Pl

    // False-asymmetry lambda0(phi):
    double a1, b1, a2, b2;
    GetFalseAsymmetryCoeffs(phi, a1, b1, a2, b2);

    const double lambda0 =
        a1 * std::cos(phi) + b1 * std::sin(phi)
      + a2 * std::cos(2.0*phi) + b2 * std::sin(2.0*phi);
      // + ... higher order terms, later
      
    // Store event
    EventCoeff ev;
    ev.lambda0 = lambda0;
    ev.lam[0]  = lambda_eff[0];
    ev.lam[1]  = lambda_eff[1];
    ev.lam[2]  = lambda_eff[2];
    ev.weight  = 1.0;

    g_events.push_back(ev);
  }

  std::cout << "[BuildEventList] Filled " << g_events.size()
            << " events after cuts\n";
}

// ----------------------------------------------------------------------------
// 4) Minuit FCN: negative log-likelihood for (Pt, Pn, Pl)
// ----------------------------------------------------------------------------
void pol_ML_fcn(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval, Double_t* par, Int_t /*iflag*/)
{
  const double Pt = par[0];
  const double Pn = par[1];
  const double Pl = par[2];

  double nll = 0.0;

  for (const auto& ev : g_events) {
    const double beta = ev.lam[0]*Pt + ev.lam[1]*Pn + ev.lam[2]*Pl;
    const double arg  = 1.0 + ev.lambda0 + beta;

    if (arg <= 0.0) {
      // Penalize unphysical region where pdf <= 0
      nll += ev.weight * (1e6 + 1e3 * std::fabs(arg));
      continue;
    }

    nll -= ev.weight * std::log(arg);
  }

  fval = nll;
}

// ----------------------------------------------------------------------------
// 5) Driver: load file/tree, build events, run MLE, print results
// ----------------------------------------------------------------------------
void pol_transfer_MLE(const char* input_file = "hist/hist_genrp_1000.root",
                      const char* tree_name  = "Tout")
{
  TFile* f = TFile::Open(input_file, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[pol_transfer_ML] ERROR: cannot open " << input_file << "\n";
    return;
  }

  TTree* T = dynamic_cast<TTree*>(f->Get(tree_name));
  if (!T) {
    std::cerr << "[pol_transfer_ML] ERROR: cannot find tree " << tree_name
              << " in file.\n";
    return;
  }

  BuildEventList(T);

  if (g_events.empty()) {
    std::cerr << "[pol_transfer_ML] No events selected. Nothing to fit.\n";
    return;
  }

  TMinuit minuit(3);
  minuit.SetFCN(pol_ML_fcn);
  minuit.SetPrintLevel(1);

  // Initial guesses
  const double Pt0 = 0.0;
  const double Pn0 = 0.0;
  const double Pl0 = 0.0;

  minuit.DefineParameter(0, "Pt", Pt0, 0.01, -1.0, 1.0);
  minuit.DefineParameter(1, "Pn", Pn0, 0.01, -1.0, 1.0);
  minuit.DefineParameter(2, "Pl", Pl0, 0.01, -1.0, 1.0);

  minuit.Migrad();
  minuit.mnimpr();

  double Pt_best, Pn_best, Pl_best;
  double dPt, dPn, dPl;

  minuit.GetParameter(0, Pt_best, dPt);
  minuit.GetParameter(1, Pn_best, dPn);
  minuit.GetParameter(2, Pl_best, dPl);

  std::cout << "\n======== Maximum-Likelihood Results ========\n";
  std::cout << "Pt = " << Pt_best << " +/- " << dPt << "\n";
  std::cout << "Pn = " << Pn_best << " +/- " << dPn << "\n";
  std::cout << "Pl = " << Pl_best << " +/- " << dPl << "\n";

  // Covariance matrix
  TMatrixD cov(3,3);
  minuit.mnemat(cov.GetMatrixArray(), 3);

  std::cout << "\nCovariance matrix (Pt, Pn, Pl):\n";
  cov.Print();

  // Ratio Pt/Pl with proper covariance propagation (INCLUDING correlation)
  if (Pl_best != 0.0) {
    const double ratio = Pt_best / Pl_best;

    // var(r) = (dr/dPt)^2 var(Pt) + (dr/dPl)^2 var(Pl) + 2 (dr/dPt)(dr/dPl) cov(Pt,Pl)
    // dr/dPt = 1/Pl, dr/dPl = -Pt/Pl^2 = -r/Pl
    const double d_rdPt = 1.0 / Pl_best;
    const double d_rdPl = -Pt_best / (Pl_best*Pl_best);

    const double var_r =
        d_rdPt*d_rdPt * cov(0,0)
      + d_rdPl*d_rdPl * cov(2,2)
      + 2.0*d_rdPt*d_rdPl * cov(0,2);

    const double dr = (var_r > 0.0) ? std::sqrt(var_r) : 0.0;

    std::cout << "\nPt/Pl = " << ratio << " +/- " << dr << "  (with cov)\n";
  } else {
    std::cout << "\nPl ~ 0, cannot form Pt/Pl.\n";
  }

  std::cout << "===========================================\n";
}
