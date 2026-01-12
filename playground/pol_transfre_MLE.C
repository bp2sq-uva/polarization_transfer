// pol_transfer_ML.C
//
// Maximum-likelihood extraction of 3 polarization components
//   P_t, P_n, P_l
// following the structure of Puckett's thesis.
//
// You MUST fill in the marked TODO sections with how you get:
//   - beam helicity ε_i
//   - beam polarization h
//   - analyzing power A_y(p, θ)
//   - spin-transport matrix elements S_xj, S_yj
//   - FPP azimuth φ_i
//
// This is written as a ROOT macro.

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMinuit.h"
#include "TMatrixD.h"

#include <vector>
#include <cmath>
#include <iostream>

// ----------------------------------------------------------------------------
// 1. Event coefficients stored for the MLE
// ----------------------------------------------------------------------------

struct EventCoeff {
  double alpha;      // false-asymmetry combination: a1 + a2 cos2φ + b2 sin2φ
  double lam[3];     // tilde-lambda: h * ε * λ_j (j = 0:t, 1:n, 2:l)
  double weight;     // event weight (usually 1.0)
};

// Global container of events used by the FCN:
std::vector<EventCoeff> g_events;

// (optional) global beam polarization if same for all events:
double g_beam_polarization = 0.85;  // <-- TODO: set properly

// ----------------------------------------------------------------------------
// 2. Physics helper stubs that YOU must implement / adapt
// ----------------------------------------------------------------------------

// Example: get beam helicity ε_i = ±1 for event i
// You will likely have a helicity branch in genrp_tree or your replay tree.
int GetHelicityForEvent(/* your branch accessors here */)
{
  // TODO: replace with your actual helicity logic
  // Example placeholder:
  // return (helicity_branch_value > 0 ? +1 : -1);
  return +1;
}

// Example: analyzing power parameterization A_y(p, θ)
// For now, you can return 1.0 to debug the framework, then plug in the real Ay.
double GetAnalyzingPower(double p, double th)
{
  // TODO: implement your Ay(p, θ) (from fit, parameterization or lookup table)
  // For testing:
  return 1.0;
}

// Example: compute FPP scattering azimuth φ for this event
// In practice, you get this from FPP track coordinates in the usual FPP frame.
double GetFPPphi(/* your branch accessors here */)
{
  // TODO: obtain φ from your FPP track: atan2(y', x') or whatever convention you use
  // Return radians.
  return 0.0;
}

// Example: compute spin-transport matrix elements from target → focal plane
// We only need S_xj and S_yj (j = t, n, l) to build λ_j.
// You likely have a COSY-based function already; call that here.
void ComputeSpinTransportForEvent(
    /* your event-level kinematics / target coords */,
    double Sx[3],  // out: S_xt, S_xn, S_xl
    double Sy[3]   // out: S_yt, S_yn, S_yl
)
{
  // TODO: call your COSY spin-transport code here, something like:
  //
  //   SpinMatrix S = GetSpinMatrixFromCOSY(x_tar, y_tar, x'_tar, y'_tar, δ, ...);
  //   Sx[0] = S.Sxt;  Sx[1] = S.Sxn;  Sx[2] = S.Sxl;
  //   Sy[0] = S.Syt;  Sy[1] = S.Syn;  Sy[2] = S.Syl;
  //
  // For now, identity approx (no precession) just for testing:
  Sx[0] = 0.0; Sx[1] = 0.0; Sx[2] = 0.0;
  Sy[0] = 1.0; Sy[1] = 0.0; Sy[2] = 0.0;
}

// False-asymmetry Fourier coefficients a1, a2, b2
// In Puckett's analysis, these are usually obtained from helicity-summed φ
// distributions and then fixed. Here you can either:
//  - pass them in as constants, or
//  - read them from a config, or
//  - set them to zero initially.
void GetFalseAsymmetryCoeffs(double phi, double& a1, double& a2, double& b2)
{
  // TODO: implement your actual false-asymmetry model (from previous fits).
  // For now, start with zero to test the framework:
  a1 = 0.0;
  a2 = 0.0;
  b2 = 0.0;
}

// ----------------------------------------------------------------------------
// 3. Build the event list g_events from a TTree (or TChain)
// ----------------------------------------------------------------------------
//
// This example assumes the TTree is called "T" and lives in one file.
// Replace the branch access comments with YOUR actual branches.
// ----------------------------------------------------------------------------

void BuildEventList(TTree* T)
{
  g_events.clear();

  // ------------------------------------------------------------------------
  // TODO: Set branch addresses for the stuff you actually have:
  //   - FPP kinematics (to get p, θ, φ)
  //   - helicity
  //   - any variables needed for COSY spin matrix
  // ------------------------------------------------------------------------

  // Example placeholders (you MUST adapt these to your actual branches):
  double p_fpp = 0.0;        // secondary scattered proton momentum
  double th_fpp = 0.0;       // secondary polar angle at FPP
  // double x_fpp, y_fpp, ... etc.

  // Branch setup examples:
  // T->SetBranchAddress("fpp.p", &p_fpp);
  // T->SetBranchAddress("fpp.th", &th_fpp);
  // T->SetBranchAddress("fpp.x", &x_fpp);
  // T->SetBranchAddress("fpp.y", &y_fpp);
  // T->SetBranchAddress("helicity", &helicity_branch_value);
  //
  // Also set whatever you need for spin transport:
  //   x_tar, y_tar, x'_tar, y'_tar, δ, ...

  Long64_t nentries = T->GetEntries();
  std::cout << "[BuildEventList] Processing " << nentries << " events\n";

  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);

    // Apply your usual elastic + quality cuts:
    // if (!PassElasticCut()) continue;
    // if (!PassFPPTrackCut()) continue;

    // Get helicity ε_i:
    int eps = GetHelicityForEvent(/* use your branch values here */);

    // Beam polarization h (if event-by-event you could pass that too)
    double h = g_beam_polarization;

    // FPP kinematics, analyzing power, φ:
    double phi = GetFPPphi(/*use your FPP branches here*/);
    double Ay  = GetAnalyzingPower(p_fpp, th_fpp);

    // Spin transport S_xj, S_yj:
    double Sx[3], Sy[3];
    ComputeSpinTransportForEvent(
        /* your event kinematics here */, Sx, Sy);

    // Build λ_j for this event:
    // λ_j = A_y ( S_yj cos φ - S_xj sin φ )
    double lambda_raw[3];
    double c = std::cos(phi);
    double s = std::sin(phi);
    for (int j = 0; j < 3; j++) {
      lambda_raw[j] = Ay * ( Sy[j]*c - Sx[j]*s );
    }

    // Fold in helicity and beam polarization: tilde-λ_j = h * ε * λ_j
    double tilde_lam[3];
    double sign = h * eps;
    for (int j = 0; j < 3; j++) {
      tilde_lam[j] = sign * lambda_raw[j];
    }

    // False-asymmetry combination α_i = a1 + a2 cos2φ + b2 sin2φ:
    double a1, a2, b2;
    GetFalseAsymmetryCoeffs(phi, a1, a2, b2);
    double alpha = a1 + a2 * std::cos(2.0*phi) + b2 * std::sin(2.0*phi);

    // Build the EventCoeff struct:
    EventCoeff ev;
    ev.alpha   = alpha;
    ev.lam[0]  = tilde_lam[0]; // t
    ev.lam[1]  = tilde_lam[1]; // n
    ev.lam[2]  = tilde_lam[2]; // l
    ev.weight  = 1.0;          // change if you have special weights

    g_events.push_back(ev);
  }

  std::cout << "[BuildEventList] Filled " << g_events.size()
            << " events after cuts\n";
}

// ----------------------------------------------------------------------------
// 4. Minuit FCN: negative log-likelihood for (Pt, Pn, Pl)
// ----------------------------------------------------------------------------
//
// Parameters:
//    par[0] = Pt
//    par[1] = Pn
//    par[2] = Pl
// ----------------------------------------------------------------------------

void pol_ML_fcn(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag)
{
  double Pt = par[0];
  double Pn = par[1];
  double Pl = par[2];

  double nll = 0.0;

  for (const auto& ev : g_events) {
    // β_i = Σ_j lam_ij * P_j
    double beta = ev.lam[0]*Pt + ev.lam[1]*Pn + ev.lam[2]*Pl;
    double arg  = 1.0 + ev.alpha + beta;

    if (arg <= 0.0) {
      // Penalize unphysical region where probability ~ arg <= 0
      nll += ev.weight * (1e6 + 1e3 * std::fabs(arg));
      continue;
    }

    nll -= ev.weight * std::log(arg);
  }

  fval = nll;
}

// ----------------------------------------------------------------------------
// 5. Driver: load file/tree, build events, run MLE, print results
// ----------------------------------------------------------------------------

void pol_transfer_MLE(const char* input_file = "genrp_data.root",
                     const char* tree_name  = "T")
{
  // Open file and get TTree
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

  // 1) Build event list for MLE
  BuildEventList(T);

  if (g_events.empty()) {
    std::cerr << "[pol_transfer_ML] No events selected. Nothing to fit.\n";
    return;
  }

  // 2) Set up Minuit with 3 parameters: Pt, Pn, Pl
  TMinuit minuit(3);
  minuit.SetFCN(pol_ML_fcn);
  minuit.SetPrintLevel(1);  // 0: quiet, 1: normal

  // Initial guesses (can use Born / super-ratio / theory)
  double Pt0 = 0.0;
  double Pn0 = 0.0;
  double Pl0 = 0.0;

  // Define parameters: name, initial value, step, low, high
  minuit.DefineParameter(0, "Pt", Pt0, 0.01, -1.0, 1.0);
  minuit.DefineParameter(1, "Pn", Pn0, 0.01, -1.0, 1.0);
  minuit.DefineParameter(2, "Pl", Pl0, 0.01, -1.0, 1.0);

  // If you want to fix Pn = 0 (like ep), uncomment:
  // minuit.FixParameter(1);

  int ierr = 0;
  minuit.Migrad();   // main minimization
  minuit.mnimpr();   // optional improvement

  // 3) Extract best-fit parameters and uncertainties
  double Pt_best, Pn_best, Pl_best;
  double dPt, dPn, dPl;

  minuit.GetParameter(0, Pt_best, dPt);
  minuit.GetParameter(1, Pn_best, dPn);
  minuit.GetParameter(2, Pl_best, dPl);

  std::cout << "\n======== Maximum-Likelihood Results ========\n";
  std::cout << "Pt = " << Pt_best << " +/- " << dPt << "\n";
  std::cout << "Pn = " << Pn_best << " +/- " << dPn << "\n";
  std::cout << "Pl = " << Pl_best << " +/- " << dPl << "\n";

  // 4) Covariance matrix (and correlations if you want)
  TMatrixD cov(3,3);
  minuit.mnemat(cov.GetMatrixArray(), 3);

  std::cout << "\nCovariance matrix (Pt, Pn, Pl):\n";
  cov.Print();

  // 5) Example: ratio Pt/Pl with error propagation (ignoring correlations)
  if (Pl_best != 0.0) {
    double ratio = Pt_best / Pl_best;
    // naive error propagation (no correlations):
    double dr = std::fabs(ratio) * std::sqrt( (dPt/Pt_best)*(dPt/Pt_best)
                                             + (dPl/Pl_best)*(dPl/Pl_best) );
    std::cout << "\nPt/Pl = " << ratio << " +/- " << dr << "  (no corr)\n";
  } else {
    std::cout << "\nPl ~ 0, cannot form Pt/Pl.\n";
  }

  std::cout << "===========================================\n";
}
