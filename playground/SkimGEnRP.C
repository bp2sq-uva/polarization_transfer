#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <fstream>

using namespace std;


//-----------------------------------------------------------------------------------------------------------------------------

void SkimGEnRP( Int_t run_no = 946 ) {

  TChain* C = new TChain("T");
//   C->Add(Form("/volatile/halla/sbs/bhasitha/Replays/Replays_Bhasitha/ToF_Calib_val_LD2/rootfiles/e1217004_fullreplay_%d_*.root", run_no));
  C->Add("/volatile/halla/sbs/bhasitha/Replays/Replays_Bhasitha/ToF_Calib_val_LD2/rootfiles/e1217004_fullreplay_*.root");

  TTreeReader Tl(C);

  TTreeReaderValue<double> scalhel_hel(Tl, "scalhel.hel");
  double scalhel_helr;

  TTreeReaderValue<double> bb_tr_n(Tl, "bb.tr.n");
  TTreeReaderValue<double> sbs_gemCeF_track_ntrack(Tl, "sbs.gemCeF.track.ntrack");  
  TTreeReaderValue<double> sbs_gemCeR_track_ntrack(Tl, "sbs.gemCeR.track.ntrack");  
  double bb_tr_nr;
  double sbs_gemCeF_track_ntrackr;
  double sbs_gemCeR_track_ntrackr;

  TTreeReaderValue<double> bb_ps_e(Tl, "bb.ps.e");
  TTreeReaderValue<double> bb_sh_e(Tl, "bb.sh.e");
  TTreeReaderValue<double> sbs_hcal_x(Tl, "sbs.hcal.x");
  TTreeReaderValue<double> sbs_hcal_y(Tl, "sbs.hcal.y");
  TTreeReaderValue<double> sbs_hcal_atimeblk(Tl, "sbs.hcal.atimeblk");
  TTreeReaderValue<double> bb_sh_atimeblk(Tl, "bb.sh.atimeblk");
  double bb_ps_er;
  double bb_sh_er;
  double sbs_hcal_xr;
  double sbs_hcal_yr;
  double sbs_hcal_atimeblkr;
  double bb_sh_atimeblkr;

  TTreeReaderValue<double> e_kine_W2(Tl, "e.kine.W2");
  double e_kine_W2r;

  TTreeReaderArray<double> bb_tr_px(Tl, "bb.tr.px");
  TTreeReaderArray<double> bb_tr_py(Tl, "bb.tr.py");
  TTreeReaderArray<double> bb_tr_pz(Tl, "bb.tr.pz");
  TTreeReaderArray<double> bb_tr_vx(Tl, "bb.tr.vx");
  TTreeReaderArray<double> bb_tr_vy(Tl, "bb.tr.vy");
  TTreeReaderArray<double> bb_tr_vz(Tl, "bb.tr.vz");
  TTreeReaderArray<double> bb_tr_p(Tl, "bb.tr.p");
  TTreeReaderArray<double> bb_tr_x(Tl, "bb.tr.x");
  TTreeReaderArray<double> bb_tr_y(Tl, "bb.tr.y");
  TTreeReaderArray<double> bb_tr_th(Tl, "bb.tr.th");
  TTreeReaderArray<double> bb_tr_ph(Tl, "bb.tr.ph");
  double bb_tr_pxr;
  double bb_tr_pyr;
  double bb_tr_pzr;
  double bb_tr_vxr;
  double bb_tr_vyr;
  double bb_tr_vzr;
  double bb_tr_pr;
  double bb_tr_xr;
  double bb_tr_yr;
  double bb_tr_thr;
  double bb_tr_phr;

  TTreeReaderArray<double> bb_tr_tg_th(Tl, "bb.tr.tg_th");
  TTreeReaderArray<double> bb_tr_tg_ph(Tl, "bb.tr.tg_ph");
  TTreeReaderArray<double> bb_tr_tg_y(Tl, "bb.tr.tg_y");
  double bb_tr_tg_thr;
  double bb_tr_tg_phr;
  double bb_tr_tg_yr;

  TTreeReaderArray<double> sbs_gemCeR_track_x(Tl, "sbs.gemCeR.track.x");
  TTreeReaderArray<double> sbs_gemCeR_track_y(Tl, "sbs.gemCeR.track.y");
  TTreeReaderArray<double> sbs_gemCeR_track_xp(Tl, "sbs.gemCeR.track.xp");
  TTreeReaderArray<double> sbs_gemCeR_track_yp(Tl, "sbs.gemCeR.track.yp");
  TTreeReaderArray<double> sbs_gemCeR_track_sclose(Tl, "sbs.gemCeR.track.sclose");
  TTreeReaderArray<double> sbs_gemCeR_track_zclose(Tl, "sbs.gemCeR.track.zclose");
  TTreeReaderArray<double> sbs_gemCeR_track_theta(Tl, "sbs.gemCeR.track.theta");
  TTreeReaderArray<double> sbs_gemCeR_track_phi(Tl, "sbs.gemCeR.track.phi");
  double sbs_gemCeR_track_xr;
  double sbs_gemCeR_track_yr; 
  double sbs_gemCeR_track_xpr;
  double sbs_gemCeR_track_ypr;
  double sbs_gemCeR_track_scloser;
  double sbs_gemCeR_track_zcloser;
  double sbs_gemCeR_track_thetar;
  double sbs_gemCeR_track_phir;

  //-----------------------------------------------------------------------------------------------------------------------------

  const Int_t nphibins = 36;

  TFile *outfile = new TFile(Form("hist/skimh/skim_genrp_%i.root",run_no),"RECREATE");

  TTree* outtree = new TTree("TSkim","Output Tree"); 
  outtree->SetAutoSave(); 

  outtree->Branch("scalhel_hel", &scalhel_helr, "scalhel_hel/D");   

  outtree->Branch("bb_tr_n", &bb_tr_nr, "bb_tr_n/D");   
  outtree->Branch("sbs_gemCeF_track_ntrack", &sbs_gemCeF_track_ntrackr, "sbs_gemCeF_track_ntrack/D");   
  outtree->Branch("sbs_gemCeR_track_ntrack", &sbs_gemCeR_track_ntrackr, "sbs_gemCeR_track_ntrack/D");   

  outtree->Branch("bb_ps_e", &bb_ps_er, "bb_ps_e/D");   
  outtree->Branch("bb_sh_e", &bb_sh_er, "bb_sh_e/D");   
  outtree->Branch("sbs_hcal_x", &sbs_hcal_xr, "sbs_hcal_x/D");   
  outtree->Branch("sbs_hcal_y", &sbs_hcal_yr, "sbs_hcal_y/D");   
  outtree->Branch("sbs_hcal_atimeblk", &sbs_hcal_atimeblkr, "sbs_hcal_atimeblk/D");   
  outtree->Branch("bb_sh_atimeblk", &bb_sh_atimeblkr, "bb_sh_atimeblk/D");   

  outtree->Branch("e_kine_W2", &e_kine_W2r, "e_kine_W2/D");   
  outtree->Branch("bb_tr_px", &bb_tr_pxr, "bb_tr_px/D");   
  outtree->Branch("bb_tr_py", &bb_tr_pyr, "bb_tr_py/D");   
  outtree->Branch("bb_tr_pz", &bb_tr_pzr, "bb_tr_pz/D");   
  outtree->Branch("bb_tr_vx", &bb_tr_vxr, "bb_tr_vx/D");   
  outtree->Branch("bb_tr_vy", &bb_tr_vyr, "bb_tr_vy/D");   
  outtree->Branch("bb_tr_vz", &bb_tr_vzr, "bb_tr_vz/D");   
  outtree->Branch("bb_tr_p", &bb_tr_pr, "bb_tr_p/D");   
  outtree->Branch("bb_tr_x", &bb_tr_xr, "bb_tr_x/D");   
  outtree->Branch("bb_tr_y", &bb_tr_yr, "bb_tr_y/D");   
  outtree->Branch("bb_tr_th", &bb_tr_thr, "bb_tr_th/D");   
  outtree->Branch("bb_tr_ph", &bb_tr_phr, "bb_tr_ph/D");   

  outtree->Branch("bb_tr_tg_th", &bb_tr_tg_thr, "bb_tr_tg_th/D");   
  outtree->Branch("bb_tr_tg_ph", &bb_tr_tg_phr, "bb_tr_tg_ph/D");   
  outtree->Branch("bb_tr_tg_y", &bb_tr_tg_yr, "bb_tr_tg_y/D");   

  outtree->Branch("sbs_gemCeR_track_x", &sbs_gemCeR_track_xr, "sbs_gemCeR_track_x/D");   
  outtree->Branch("sbs_gemCeR_track_y", &sbs_gemCeR_track_yr, "sbs_gemCeR_track_y/D");   
  outtree->Branch("sbs_gemCeR_track_xp", &sbs_gemCeR_track_xpr, "sbs_gemCeR_track_xp/D");   
  outtree->Branch("sbs_gemCeR_track_yp", &sbs_gemCeR_track_ypr, "sbs_gemCeR_track_yp/D");   
  outtree->Branch("sbs_gemCeR_track_sclose", &sbs_gemCeR_track_scloser, "sbs_gemCeR_track_sclose/D");   
  outtree->Branch("sbs_gemCeR_track_zclose", &sbs_gemCeR_track_zcloser, "sbs_gemCeR_track_zclose/D");   
  outtree->Branch("sbs_gemCeR_track_theta", &sbs_gemCeR_track_thetar, "sbs_gemCeR_track_theta/D");   
  outtree->Branch("sbs_gemCeR_track_phi", &sbs_gemCeR_track_phir, "sbs_gemCeR_track_phi/D");   

  //-----------------------------------------------------------------------------------------------------------------------------

  Long64_t ev = 0;

  while( Tl.Next() ) {

    scalhel_helr = *scalhel_hel;

    bb_tr_nr = *bb_tr_n;
    sbs_gemCeF_track_ntrackr = *sbs_gemCeF_track_ntrack;
    sbs_gemCeR_track_ntrackr = *sbs_gemCeR_track_ntrack;

    bb_ps_er = *bb_ps_e;
    bb_sh_er = *bb_sh_e;
    sbs_hcal_xr = *sbs_hcal_x;
    sbs_hcal_yr = *sbs_hcal_y;
    sbs_hcal_atimeblkr = *sbs_hcal_atimeblk;
    bb_sh_atimeblkr = *bb_sh_atimeblk;

    e_kine_W2r = *e_kine_W2;

    bb_tr_pxr = bb_tr_px[0];
    bb_tr_pyr = bb_tr_py[0];
    bb_tr_pzr = bb_tr_pz[0];
    bb_tr_pzr = bb_tr_pz[0];
    bb_tr_vxr = bb_tr_vx[0];
    bb_tr_vyr = bb_tr_vy[0];
    bb_tr_vzr = bb_tr_vz[0];
    bb_tr_pr  = bb_tr_p[0];
    bb_tr_xr  = bb_tr_x[0];
    bb_tr_yr  = bb_tr_y[0];
    bb_tr_thr = bb_tr_th[0];
    bb_tr_phr = bb_tr_ph[0];

    bb_tr_tg_thr = bb_tr_tg_th[0];
    bb_tr_tg_phr = bb_tr_tg_ph[0];
    bb_tr_tg_yr  = bb_tr_tg_y[0];

    sbs_gemCeR_track_xr = sbs_gemCeR_track_x[0];
    sbs_gemCeR_track_yr = sbs_gemCeR_track_y[0];
    sbs_gemCeR_track_xpr = sbs_gemCeR_track_xp[0];
    sbs_gemCeR_track_ypr = sbs_gemCeR_track_yp[0];
    sbs_gemCeR_track_scloser = sbs_gemCeR_track_sclose[0];
    sbs_gemCeR_track_zcloser = sbs_gemCeR_track_zclose[0];
    sbs_gemCeR_track_thetar = sbs_gemCeR_track_theta[0];
    sbs_gemCeR_track_phir = sbs_gemCeR_track_phi[0];
    
    outtree->Fill();

    ev++;
  }
  
  outfile->Write();
  outfile->Close();

}

//-----------------------------------------------------------------------------------------------------------------------------