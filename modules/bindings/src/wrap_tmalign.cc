//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

#include "tmalign/TMalign.h" // include for the external TMalign

#include <ost/mol/atom_view.hh>
#include <ost/message.hh>
#include <ost/bindings/wrap_tmalign.hh>

namespace ost{ namespace bindings{

TMAlignResult WrappedTMAlign(const geom::Vec3List& pos_one, 
                             const geom::Vec3List& pos_two, 
                             const ost::seq::SequenceHandle& seq1,
                             const ost::seq::SequenceHandle& seq2,
                             bool fast) {

  int xlen = pos_one.size();
  int ylen = pos_two.size();  

  if(xlen <= 5 || ylen <= 5) {
    throw ost::Error("Input sequence too short!");
  }

  if(xlen != seq1.GetLength() || ylen != seq2.GetLength()) {
    throw ost::Error("Positions and Sequence must have consistent size "
                     "to run TMAlign!");
  }

  // squeeze input into right format
  char* seqx = new char[xlen+1];
  char* seqy = new char[ylen+1];
  seqx[xlen] = '\0';
  seqy[ylen] = '\0';
  char* secx = new char[xlen];
  char* secy = new char[ylen];

  // use TMalign functionality to generate position arrays
  double** xa;
  double** ya;
  NewArray(&xa, xlen, 3);
  NewArray(&ya, ylen, 3);

  for(int i = 0; i < xlen; ++i) {
    xa[i][0] = pos_one[i][0];
    xa[i][1] = pos_one[i][1];
    xa[i][2] = pos_one[i][2];
    seqx[i] = seq1[i];
  }

  for(int i = 0; i < ylen; ++i) {
    ya[i][0] = pos_two[i][0];
    ya[i][1] = pos_two[i][1];
    ya[i][2] = pos_two[i][2];
    seqy[i] = seq2[i];
  }

  make_sec(xa, xlen, secx);
  make_sec(ya, ylen, secy);

  // these variables are chosen such that running TMalign_main is the same as 
  // you would call the executable without any additional parameters
  double Lnorm_ass = -1.0;
  double d0_scale = -1.0;
  bool u_opt = false;
  bool d_opt = false;
  int a_opt = 0; 
  std::vector<String> sequence; 
  bool i_opt = false;
  double TMcut = -1; 

  // following variables are copied from the TMAlign source code
  double t0[3], u0[3][3];
  double TM1, TM2;
  double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
  double d0_0, TM_0;
  double d0A, d0B, d0u, d0a;
  double d0_out=5.0;
  String seqM, seqxA, seqyA;// for output alignment
  double rmsd0 = 0.0;
  int L_ali;                // Aligned length in standard_TMscore
  double Liden=0;
  double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
  int n_ali=0;
  int n_ali8=0;

  TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0, TM1, TM2, TM3, TM4, TM5,
               d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
               rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen, 
               sequence, Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt, 
               fast, 0, TMcut);

  // cleanup
  DeleteArray(&xa, xlen);
  DeleteArray(&ya, ylen);
  delete [] seqx;
  delete [] seqy;
  delete [] secx;
  delete [] secy;

  // collect results and return
  TMAlignResult res;
  res.tm_score = TM1;
  res.rmsd = rmsd0;
  res.aligned_length = n_ali8;
  res.transform = geom::Mat4(u0[0][0], u0[0][1], u0[0][2], t0[0],
                             u0[1][0], u0[1][1], u0[1][2], t0[1],
                             u0[2][0], u0[2][1], u0[2][2], t0[2],
                             0.0, 0.0, 0.0, 1.0);
  res.alignment = ost::seq::CreateAlignment();
  ost::seq::SequenceHandle aligned_seq1 = ost::seq::CreateSequence("seq1", seqxA);
  ost::seq::SequenceHandle aligned_seq2 = ost::seq::CreateSequence("seq2", seqyA);
  res.alignment.AddSequence(aligned_seq1);
  res.alignment.AddSequence(aligned_seq2);
  return res;
}


TMAlignResult WrappedTMAlign(const ost::mol::ChainView& chain1,
                             const ost::mol::ChainView& chain2,
                             bool fast) {

  geom::Vec3List pos1;
  geom::Vec3List pos2;
  std::vector<char> s1;
  std::vector<char> s2;

  ost::mol::ResidueViewList res_list_1 = chain1.GetResidueList();
  ost::mol::ResidueViewList res_list_2 = chain2.GetResidueList();

  for(ost::mol::ResidueViewList::iterator it = res_list_1.begin(); 
      it != res_list_1.end(); ++it) {
    if(!it->IsPeptideLinking()) {
      continue;
    }
    ost::mol::AtomView ca = it->FindAtom("CA");
    if(!ca.IsValid()) {
      continue;
    }
    char olc = it->GetOneLetterCode();
    if(olc == '?') {
      continue;
    }
    pos1.push_back(ca.GetPos());
    s1.push_back(olc);
  }

  for(ost::mol::ResidueViewList::iterator it = res_list_2.begin(); 
      it != res_list_2.end(); ++it) {
    if(!it->IsPeptideLinking()) {
      continue;
    }
    ost::mol::AtomView ca = it->FindAtom("CA");
    if(!ca.IsValid()) {
      continue;
    }
    char olc = it->GetOneLetterCode();
    if(olc == '?') {
      continue;
    }
    pos2.push_back(ca.GetPos());
    s2.push_back(olc);
  }

  String str_s1(s1.begin(), s1.end());
  String str_s2(s2.begin(), s2.end());

  ost::seq::SequenceHandle seq_s1 = ost::seq::CreateSequence("one", str_s1);
  ost::seq::SequenceHandle seq_s2 = ost::seq::CreateSequence("two", str_s2);

  return WrappedTMAlign(pos1, pos2, seq_s1, seq_s2, fast);
}

}} //ns
