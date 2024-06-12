//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2024 by the OpenStructure authors
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

#include <ost/bindings/wrap_parasail.hh>

#if OST_PARASAIL_ENABLED

#include <parasail.h> // include for external parasail

typedef parasail_result_t* (*parasail_aln_fn)(const char*, int, const char*,
                                              int, int, int,
                                              const parasail_matrix_t *);

namespace{

ost::seq::AlignmentHandle Align(const ost::seq::SequenceHandle& s1,
                                const ost::seq::SequenceHandle& s2,
                                int gap_open_penalty, int gap_extension_penalty,
                                const parasail_matrix_t* matrix,
                                parasail_aln_fn aln_fn,
                                bool set_offset) {

    parasail_result_t *r = NULL;
    parasail_traceback_t *traceback_r = NULL;
    
    String s1_str = s1.GetGaplessString();
    String s2_str = s2.GetGaplessString();

    const char* s1_c_str = s1_str.c_str();
    int s1_len = s1_str.size();
    const char* s2_c_str = s2_str.c_str();
    int s2_len = s2_str.size();
    
    r = aln_fn(s1_c_str, s1_len, s2_c_str, s2_len, 11, 1, matrix);
    traceback_r = parasail_result_get_traceback(r, s1_c_str, s1_len,
                                                s2_c_str, s2_len, matrix,
                                                '|', ':', '.');
    
    ost::seq::AlignmentHandle aln = ost::seq::CreateAlignment();
    ost::seq::SequenceHandle aln_s1 =
    ost::seq::CreateSequence(s1.GetName(), String(traceback_r->query));
    
    ost::seq::SequenceHandle aln_s2 =
    ost::seq::CreateSequence(s2.GetName(), String(traceback_r->ref));

    if(set_offset) {

      /*
      Tried the code below to set offsets but observed wrongly reported
      offset which can be triggered from Python:

      from ost import io
      from ost import seq
      from ost import bindings

      test_seq = seq.CreateSequence("s", "MMDKHKYRVEIQQMMFVSGEINDPPVETTSLIEDIVRGQVIEILLQSNKTAHLRGSRSIL"
                                         "PEDVIFLIRHDKAKVNRLRTYLSWKDLRKNAKDQDASAGVASGTGNPGAGGEDDLKKAGG"
                                         "GEKDEKDGGNMMKVKKSQIKLPWELQFMFNEHPLENNDDNDDMDEDEREANIVTLKRLKM"
                                         "ADDRTRNMTKEEYVHWSDCRQASFTFRKNKRFKDWSGISQLTEGKPHDDVIDILGFLTFE"
                                         "IVCSLTETALKIKQREQVLQTQKDKSQQSSQDNTNFEFASSTLHRKKRLFDGPENVINPL"
                                         "KPRHIEEAWRVLQTIDMRHRALTNFKGGRLSSKPIIM")
                                   

      s = seq.CreateSequence("asdf", "MTKSIYIIIGYMLHDEEFFYFFFISFYTLWIVFFLLHLSFFSTLSFGIFHDFDTDVYIKVGNYIL"
                                     "HFLELSKNSNLLKNSSEMLKHFRLASLMYMYVYTQMICPSLLGIRN")


      aln = bindings.ParaLocalAlign(test_seq, s)

      print(aln.ToString())
      print(aln.GetSequence(0).offset)
      print(aln.GetSequence(1).offset)




      // Not sure if this is the most efficient way to set offsets,
      // i.e. the index of start characters which are actually
      // present in returned alignments
      parasail_cigar_t* cigar = parasail_result_get_cigar(r, s1_c_str, s1_len,
                                                          s2_c_str, s2_len,
                                                          matrix);

      aln_s1.SetOffset(cigar->beg_query);
      aln_s2.SetOffset(cigar->beg_ref);
      parasail_cigar_free(cigar);
      */

      // this is far from optimal but doesn't have the issue described above
      aln_s1.SetOffset(s1.GetString().find(aln_s1.GetGaplessString()));
      aln_s2.SetOffset(s2.GetString().find(aln_s2.GetGaplessString()));
    }
    aln.AddSequence(aln_s1);
    aln.AddSequence(aln_s2);
    
    parasail_result_free(r);
    parasail_traceback_free(traceback_r);
    
    return aln;
}

}

namespace ost{ namespace bindings{

ost::seq::AlignmentHandle ParaLocalAlign(const ost::seq::SequenceHandle& s1,
                                         const ost::seq::SequenceHandle& s2,
                                         int gap_open_penalty,
                                         int gap_extension_penalty,
                                         const String& matrix) {
  const parasail_matrix_t* m = parasail_matrix_lookup(matrix.c_str());
  return Align(s1, s2, gap_open_penalty, gap_extension_penalty,
               m, parasail_sw_trace_scan_sat, true); 
}

ost::seq::AlignmentHandle ParaGlobalAlign(const ost::seq::SequenceHandle& s1,
                                          const ost::seq::SequenceHandle& s2,
                                          int gap_open_penalty,
                                          int gap_extension_penalty,
                                          const String& matrix) {
  const parasail_matrix_t* m = parasail_matrix_lookup(matrix.c_str());
  return Align(s1, s2, gap_open_penalty, gap_extension_penalty,
               m, parasail_nw_trace_scan_sat, false);  
}

ost::seq::AlignmentHandle ParaSemiGlobalAlign(const ost::seq::SequenceHandle& s1,
                                              const ost::seq::SequenceHandle& s2,
                                              int gap_open_penalty,
                                              int gap_extension_penalty,
                                              const String& matrix) {
  const parasail_matrix_t* m = parasail_matrix_lookup(matrix.c_str());
  return Align(s1, s2, gap_open_penalty, gap_extension_penalty,
               m, parasail_sg_trace_scan_sat, false);   
}

bool ParasailAvailable() {
  return true;
}

}} // ns

#else

namespace ost{ namespace bindings{

ost::seq::AlignmentHandle ParaLocalAlign(const ost::seq::SequenceHandle& s1,
                                         const ost::seq::SequenceHandle& s2,
                                         int gap_open_penalty,
                                         int gap_extension_penalty,
                                         const String& matrix) {
  throw ost::Error("Must enable Parasail when compiling OpenStructure to use "
                   "ParaLocalAlign");
}

ost::seq::AlignmentHandle ParaGlobalAlign(const ost::seq::SequenceHandle& s1,
                                          const ost::seq::SequenceHandle& s2,
                                          int gap_open_penalty,
                                          int gap_extension_penalty,
                                          const String& matrix) {
  throw ost::Error("Must enable Parasail when compiling OpenStructure to use "
                   "ParaGlobalAlign");  
}

ost::seq::AlignmentHandle ParaSemiGlobalAlign(const ost::seq::SequenceHandle& s1,
                                              const ost::seq::SequenceHandle& s2,
                                              int gap_open_penalty,
                                              int gap_extension_penalty,
                                              const String& matrix) {
  throw ost::Error("Must enable Parasail when compiling OpenStructure to use "
                   "ParaSemiGlobalAlign");  
}

bool ParasailAvailable() {
  return false;
}

}} // ns

#endif
