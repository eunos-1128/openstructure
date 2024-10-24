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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/conop/conop.hh>
#include <ost/conop/heuristic.hh>
#include <ost/io/mol/mmcif_writer.hh>
#include <ost/mol/mol.hh>
#include <ost/platform.hh>
#include <ost/log.hh>

using namespace ost;
using namespace ost::conop;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

CompoundLibPtr load_lib()
{
  if (!getenv("OST_ROOT")) {
    LOG_ERROR("OST_ROOT environment variable not set. Can't load "
              "compound library without a proper OST_ROOT");
    return CompoundLibPtr();
  }
  SetPrefixPath(getenv("OST_ROOT"));
  String lib_path=GetSharedDataPath()+"/compounds.chemlib";
  CompoundLibPtr compound_lib=CompoundLib::Load(lib_path);
  return compound_lib;
}

BOOST_AUTO_TEST_CASE(mmcif_writer_force_hetatm)
{
  CompoundLibPtr lib = load_lib();
   if (!lib) { return; }
  BOOST_TEST_MESSAGE("  Running mmcif_force_hetatm tests...");
  /*
    Make sure that atoms set to HETATM are written as HETATM. There is some
    logic in place to deal with HETAM for mmcif_conform=false, check that this
    is working.
   */
  // Create small entity
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r1=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r1, "N", geom::Vec3(1, 1, 1), "N", 1.0, 0.0, true);
  edi.InsertAtom(r1, "C", geom::Vec3(4, 1, 2), "C", 1.0, 0.0, true);
  mol::ResidueHandle r2=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r2, "N", geom::Vec3(5, 2, 3), "N", 1.0, 0.0, true);
  edi.InsertAtom(r2, "C", geom::Vec3(4, 1, 1), "C", 1.0, 0.0, true);
  mol::ResidueHandle r3=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r3, "N", geom::Vec3(5, 2, 2), "N", 1.0, 0.0, true);
  edi.InsertAtom(r3, "C", geom::Vec3(4, 1, 0), "C", 1.0, 0.0, true);
  edi.SetChainType(ch, mol::CHAINTYPE_UNKNOWN);
  // make sure we have a proper polypeptide
  conop::HeuristicProcessor heu_proc;
  heu_proc.Process(ent);
  BOOST_CHECK_EQUAL(mol::InSequence(r1, r2), true);
  BOOST_CHECK_EQUAL(mol::InSequence(r2, r3), true);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, lib, false);
  std::stringstream out;
  writer.Write("test", out);

  String s=out.str();
  // make sure the entity is a polymer
  BOOST_CHECK_NE(s.find("_entity.id\n_entity.type\n1 polymer"), String::npos);
  // check all atom records to be HETATMs
  for(auto i: ch.GetAtomList()){
    BOOST_CHECK_NE(s.find("HETATM "+
                          i.GetElement()+" "+
                          i.GetName()+" "+
                          i.GetResidue().GetName()),
                   String::npos);
  }

  // check that ATOM is written, if HETATM is not set
  // Create small entity
  ent=mol::CreateEntity();
  edi=ent.EditXCS();
  ch=edi.InsertChain("A");
  r1=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r1, "N", geom::Vec3(1, 1, 1), "N");
  edi.InsertAtom(r1, "C", geom::Vec3(4, 1, 2), "C");
  r2=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r2, "N", geom::Vec3(5, 2, 3), "N");
  edi.InsertAtom(r2, "C", geom::Vec3(4, 1, 1), "C");
  r3=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r3, "N", geom::Vec3(5, 2, 2), "N");
  edi.InsertAtom(r3, "C", geom::Vec3(4, 1, 0), "C");
  edi.SetChainType(ch, mol::CHAINTYPE_UNKNOWN);
  // make sure we have a proper polypeptide
  heu_proc.Process(ent);
  BOOST_CHECK_EQUAL(mol::InSequence(r1, r2), true);
  BOOST_CHECK_EQUAL(mol::InSequence(r2, r3), true);

  // Create mmCIF stream
  writer=MMCifWriter();
  writer.SetStructure(ent, lib, false);
  out=std::stringstream();
  writer.Write("test", out);

  s=out.str();
  // make sure the entity is a polymer
  BOOST_CHECK_NE(s.find("_entity.id\n_entity.type\n1 polymer"), String::npos);
  // check all atom records to be ATOMs
  for(auto i: ch.GetAtomList()){
    BOOST_CHECK_NE(s.find("ATOM "+
                          i.GetElement()+" "+
                          i.GetName()+" "+
                          i.GetResidue().GetName()),
                   String::npos);
  }

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_writer_entity1)
{
  CompoundLibPtr lib = load_lib();
   if (!lib) { return; }
  BOOST_TEST_MESSAGE("  Running mmcif_writer_entity1 tests...");
  /*
    Make sure molecular entities in mmCIF files written by OST start counting
    at ID 1. This is not enforced by the mmCIF format definition, but common
    practice.
   */

  // Create tiny entity
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  mol::ResidueHandle r=edi.AppendResidue(ch, "GLY");
  edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5), "C");
  edi.SetChainType(ch, mol::CHAINTYPE_UNKNOWN);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, lib, false);
  std::stringstream out;
  writer.Write("test", out);

  // Check if entity starts with 1, either by reading mmCIF or "grep"
  String s=out.str();
  BOOST_CHECK_NE(s.find("_entity.id\n_entity.type\n1 non-polymer"),
                 String::npos);

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_writer_poly_vs_non_poly)
{
  CompoundLibPtr lib = load_lib();
   if (!lib) { return; }
  BOOST_TEST_MESSAGE("  Running mmcif_writer_poly_vs_non_poly tests...");
  /*
    Go for small polymers that are not polymer... the story of 2 amino acids (to
    be handled like RCSB as 2 separated chains) plus how the same thing works
    for nucleic acids.
   */

  // Polypeptide: 2aa in a chain are 2 separated entities for RCSB (check 1E8K)
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  edi.SetChainType(ch, mol::CHAINTYPE_POLY_PEPTIDE_L);
  // ALA
  mol::ResidueHandle r1=edi.AppendResidue(ch, "ALA");
  edi.InsertAtom(r1, "N",  geom::Vec3(44.987, 17.389, 12.362), "N");
  edi.InsertAtom(r1, "CA", geom::Vec3(45.936, 16.434, 12.890), "C");
  edi.InsertAtom(r1, "C",  geom::Vec3(47.196, 17.227, 13.152), "C");
  edi.InsertAtom(r1, "O",  geom::Vec3(47.506, 18.153, 12.401), "O");
  edi.InsertAtom(r1, "CB", geom::Vec3(46.244, 15.293, 11.961), "C");
  // PRO
  mol::ResidueHandle r2=edi.AppendResidue(ch, "PRO");
  edi.InsertAtom(r2, "N",   geom::Vec3(47.953, 16.910, 14.229), "N");
  edi.InsertAtom(r2, "CA",  geom::Vec3(47.673, 15.829, 15.187), "C");
  edi.InsertAtom(r2, "C",   geom::Vec3(46.564, 16.052, 16.233), "C");
  edi.InsertAtom(r2, "O",   geom::Vec3(46.059, 17.169, 16.417), "O");
  edi.InsertAtom(r2, "CB",  geom::Vec3(49.054, 15.755, 15.880), "C");
  edi.InsertAtom(r2, "CG",  geom::Vec3(49.357, 17.213, 16.030), "C");
  edi.InsertAtom(r2, "CD",  geom::Vec3(49.098, 17.714, 14.637), "C");
  edi.InsertAtom(r2, "OXT", geom::Vec3(46.144, 15.129, 16.950), "O");

  // Make sure that the two residues r1, r2 are actually connected
  conop::HeuristicProcessor heu_proc;
  heu_proc.Process(ent);
  BOOST_CHECK_EQUAL(mol::InSequence(r1, r2), true);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, lib, false);
  std::stringstream out;
  writer.Write("test", out);

  String s=out.str();
  // Check that the mmCIF output contains 2 non-polymer entities
  BOOST_CHECK_NE(
        s.find("loop_\n_entity.id\n_entity.type\n1 non-polymer\n2 non-polymer"),
        String::npos);
  BOOST_CHECK_NE(
             s.find("loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2"),
             String::npos);
  // Check that atoms are HETATMs since non-poly
  BOOST_CHECK_NE(s.find("HETATM N N ALA"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CA ALA"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C C ALA"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM O O ALA"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CB ALA"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM N N PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CA PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C C PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM O O PRO "), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CB PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CG PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM C CD PRO"), String::npos);
  BOOST_CHECK_NE(s.find("HETATM O OXT PRO"), String::npos);

  // Nucleic acid: 2 bases in a chain are 2 entities for RCSB (check 4K9A)
  ent=mol::CreateEntity();
  edi=ent.EditXCS();
  ch=edi.InsertChain("A");
  edi.SetChainType(ch, mol::CHAINTYPE_POLY_RN);
  // G
  r1=edi.AppendResidue(ch, "G");
  r1.SetChemClass(mol::ChemClass('R'));
  edi.InsertAtom(r1, "OP3", geom::Vec3(-19.992, -12.612, -22.535), "O");
  edi.InsertAtom(r1, "P",   geom::Vec3(-20.743, -13.990, -22.451), "P");
  edi.InsertAtom(r1, "OP1", geom::Vec3(-20.093, -14.853, -23.494), "O");
  edi.InsertAtom(r1, "OP2", geom::Vec3(-20.592, -14.510, -21.074), "O");
  edi.InsertAtom(r1, "O5'", geom::Vec3(-22.258, -13.779, -22.864), "O");
  edi.InsertAtom(r1, "C5'", geom::Vec3(-22.956, -12.621, -22.414), "C");
  edi.InsertAtom(r1, "C4'", geom::Vec3(-24.433, -12.831, -22.071), "C");
  edi.InsertAtom(r1, "O4'", geom::Vec3(-24.958, -11.618, -21.688), "O");
  edi.InsertAtom(r1, "C3'", geom::Vec3(-25.196, -13.187, -23.283), "C");
  edi.InsertAtom(r1, "O3'", geom::Vec3(-26.335, -13.801, -22.833), "O");
  edi.InsertAtom(r1, "C2'", geom::Vec3(-25.568, -11.908, -23.910), "C");
  edi.InsertAtom(r1, "O2'", geom::Vec3(-26.792, -12.071, -24.602), "O");
  edi.InsertAtom(r1, "C1'", geom::Vec3(-25.819, -11.112, -22.709), "C");
  edi.InsertAtom(r1, "N9",  geom::Vec3(-25.646, -9.689,  -22.974), "N");
  edi.InsertAtom(r1, "C8",  geom::Vec3(-26.518, -8.728,  -22.814), "C");
  edi.InsertAtom(r1, "N7",  geom::Vec3(-25.994, -7.595,  -23.153), "N");
  edi.InsertAtom(r1, "C5",  geom::Vec3(-24.781, -7.838,  -23.527), "C");
  edi.InsertAtom(r1, "C6",  geom::Vec3(-23.815, -7.062,  -23.942), "C");
  edi.InsertAtom(r1, "O6",  geom::Vec3(-23.953, -5.880,  -24.055), "O");
  edi.InsertAtom(r1, "N1",  geom::Vec3(-22.652, -7.575,  -24.269), "N");
  edi.InsertAtom(r1, "C2",  geom::Vec3(-22.453, -8.886,  -24.153), "C");
  edi.InsertAtom(r1, "N2",  geom::Vec3(-21.312, -9.419,  -24.474), "N");
  edi.InsertAtom(r1, "N3",  geom::Vec3(-23.422, -9.638,  -23.738), "N");
  edi.InsertAtom(r1, "C4",  geom::Vec3(-24.577, -9.122,  -23.422), "C");
  // A
  r2=edi.AppendResidue(ch, "A");
  r2.SetChemClass(mol::ChemClass('R'));
  edi.InsertAtom(r2, "P",   geom::Vec3(-27.082, -11.572, -26.017), "P");
  edi.InsertAtom(r2, "OP1", geom::Vec3(-28.487, -11.993, -26.469), "O");
  edi.InsertAtom(r2, "OP2", geom::Vec3(-27.228, -10.079, -25.961), "O");
  edi.InsertAtom(r2, "O5'", geom::Vec3(-26.066, -11.826, -27.005), "O");
  edi.InsertAtom(r2, "C5'", geom::Vec3(-25.894, -13.010, -27.689), "C");
  edi.InsertAtom(r2, "C4'", geom::Vec3(-24.500, -13.055, -28.234), "C");
  edi.InsertAtom(r2, "O4'", geom::Vec3(-24.168, -11.802, -28.703), "O");
  edi.InsertAtom(r2, "C3'", geom::Vec3(-23.540, -13.333, -27.140), "C");
  edi.InsertAtom(r2, "O3'", geom::Vec3(-22.424, -13.995, -27.663), "O");
  edi.InsertAtom(r2, "C2'", geom::Vec3(-23.101, -11.988, -26.770), "C");
  edi.InsertAtom(r2, "O2'", geom::Vec3(-21.831, -12.081, -26.257), "O");
  edi.InsertAtom(r2, "C1'", geom::Vec3(-23.037, -11.318, -28.071), "C");
  edi.InsertAtom(r2, "N9",  geom::Vec3(-23.277, -9.902,  -27.855), "N");
  edi.InsertAtom(r2, "C8",  geom::Vec3(-24.331, -9.311,  -27.370), "C");
  edi.InsertAtom(r2, "N7",  geom::Vec3(-24.073, -8.040,  -27.362), "N");
  edi.InsertAtom(r2, "C5",  geom::Vec3(-22.856, -7.821,  -27.812), "C");
  edi.InsertAtom(r2, "C6",  geom::Vec3(-22.091, -6.735,  -27.998), "C");
  edi.InsertAtom(r2, "N6",  geom::Vec3(-22.544, -5.535,  -27.711), "N");
  edi.InsertAtom(r2, "N1",  geom::Vec3(-20.885, -6.823,  -28.480), "N");
  edi.InsertAtom(r2, "C2",  geom::Vec3(-20.381, -8.034,  -28.799), "C");
  edi.InsertAtom(r2, "N3",  geom::Vec3(-21.167, -9.113,  -28.595), "N");
  edi.InsertAtom(r2, "C4",  geom::Vec3(-22.379, -8.993,  -28.116), "C");

  // Make sure that the two residues r1, r2 are actually connected
  /* In the RCSB entry 4K9A, G O2' and A P are linked, but we force a regular
     nucleotide link for the test
  */
  edi.Connect(r1.FindAtom("O3'"), r2.FindAtom("P"));
  BOOST_CHECK_EQUAL(mol::InSequence(r1, r2), true);

  // Create mmCIF stream
  writer=MMCifWriter();
  writer.SetStructure(ent, lib, false);
  out=std::stringstream();
  writer.Write("test", out);

  s=out.str();
  // Check that the mmCIF output contains 2 non-polymer entities
  BOOST_CHECK_NE(
        s.find("loop_\n_entity.id\n_entity.type\n1 non-polymer\n2 non-polymer"),
        String::npos);
  BOOST_CHECK_NE(
             s.find("loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2"),
             String::npos);
  // Check that atoms are HETATMs since non-poly
  for(auto i: ch.GetAtomList()){
    BOOST_CHECK_NE(s.find("HETATM "+
                          i.GetElement()+" "+
                          i.GetName()+" "+
                          i.GetResidue().GetName()),
                   String::npos);
  }

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(mmcif_writer_small_sugars)
{
  CompoundLibPtr lib = load_lib();
   if (!lib) { return; }
  BOOST_TEST_MESSAGE("  Running mmcif_writer_small_sugars tests...");
  /*
    While RCSB marks dipeptides and dinucleotides as non-ploymers, sugars are
    marked branched as soon as there are 2 connected.
   */
  // Branched: 2 sugars connected (check RCSB 3AXH)
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  edi.SetChainType(ch, mol::CHAINTYPE_OLIGOSACCHARIDE);
  // GLC
  mol::ResidueHandle r1=edi.AppendResidue(ch, "GLC");
  r1.SetChemClass(mol::ChemClass('Y'));
  edi.InsertAtom(r1, "C1", geom::Vec3(-17.103, -7.005, -18.605), "C");
  edi.InsertAtom(r1, "C2", geom::Vec3(-18.238, -7.769, -17.910), "C");
  edi.InsertAtom(r1, "C3", geom::Vec3(-18.607, -9.018, -18.701), "C");
  edi.InsertAtom(r1, "C4", geom::Vec3(-18.909, -8.623, -20.136), "C");
  edi.InsertAtom(r1, "C5", geom::Vec3(-17.692, -7.925, -20.738), "C");
  edi.InsertAtom(r1, "C6", geom::Vec3(-17.985, -7.413, -22.123), "C");
  edi.InsertAtom(r1, "O1", geom::Vec3(-15.956, -7.805, -18.503), "O");
  edi.InsertAtom(r1, "O2", geom::Vec3(-17.855, -8.167, -16.612), "O");
  edi.InsertAtom(r1, "O3", geom::Vec3(-19.758, -9.617, -18.139), "O");
  edi.InsertAtom(r1, "O4", geom::Vec3(-19.249, -9.772, -20.895), "O");
  edi.InsertAtom(r1, "O5", geom::Vec3(-17.381, -6.775, -19.989), "O");
  edi.InsertAtom(r1, "O6", geom::Vec3(-19.210, -6.695, -22.028), "O");
  // GLC
  mol::ResidueHandle r2=edi.AppendResidue(ch, "GLC");
  r2.SetChemClass(mol::ChemClass('Y'));
  edi.InsertAtom(r2, "C1", geom::Vec3(-20.076, -6.111, -23.424), "C");
  edi.InsertAtom(r2, "C2", geom::Vec3(-21.506, -5.777, -22.960), "C");
  edi.InsertAtom(r2, "C3", geom::Vec3(-22.102, -7.005, -22.273), "C");
  edi.InsertAtom(r2, "C4", geom::Vec3(-22.152, -8.101, -23.352), "C");
  edi.InsertAtom(r2, "C5", geom::Vec3(-20.802, -8.397, -24.036), "C");
  edi.InsertAtom(r2, "C6", geom::Vec3(-21.027, -9.040, -25.414), "C");
  edi.InsertAtom(r2, "O2", geom::Vec3(-21.398, -4.714, -22.066), "O");
  edi.InsertAtom(r2, "O3", geom::Vec3(-23.396, -6.675, -21.796), "O");
  edi.InsertAtom(r2, "O4", geom::Vec3(-22.738, -9.283, -22.830), "O");
  edi.InsertAtom(r2, "O5", geom::Vec3(-20.007, -7.236, -24.284), "O");
  edi.InsertAtom(r2, "O6", geom::Vec3(-21.528, -8.081, -26.325), "O");

  // Connect the two sugars
  edi.Connect(r1.FindAtom("O6"), r2.FindAtom("C1"));

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, lib, false);
  std::stringstream out;
  writer.Write("test", out);

  String s=out.str();
  // Check that the mmCIF output contains a branched entity
  BOOST_CHECK_NE(s.find("loop_\n_entity.id\n_entity.type\n1 branched"),
                 String::npos);
  // Check that atoms are HETATMs since non-poly
  for(auto i: ch.GetAtomList()){
    BOOST_CHECK_NE(s.find("HETATM "+
                          i.GetElement()+" "+
                          i.GetName()+" "+
                          i.GetResidue().GetName()),
                   String::npos);
  }

  // Non-poly: single sugar (check RCSB 1BDG)
  ent=mol::CreateEntity();
  edi=ent.EditXCS();
  ch=edi.InsertChain("A");
  edi.SetChainType(ch, mol::CHAINTYPE_OLIGOSACCHARIDE);
  // GLC
  r1=edi.AppendResidue(ch, "GLC");
  r1.SetChemClass(mol::ChemClass('Y'));
  edi.InsertAtom(r1, "C1", geom::Vec3(-17.103, -7.005, -18.605), "C");
  edi.InsertAtom(r1, "C2", geom::Vec3(-18.238, -7.769, -17.910), "C");
  edi.InsertAtom(r1, "C3", geom::Vec3(-18.607, -9.018, -18.701), "C");
  edi.InsertAtom(r1, "C4", geom::Vec3(-18.909, -8.623, -20.136), "C");
  edi.InsertAtom(r1, "C5", geom::Vec3(-17.692, -7.925, -20.738), "C");
  edi.InsertAtom(r1, "C6", geom::Vec3(-17.985, -7.413, -22.123), "C");
  edi.InsertAtom(r1, "O1", geom::Vec3(-15.956, -7.805, -18.503), "O");
  edi.InsertAtom(r1, "O2", geom::Vec3(-17.855, -8.167, -16.612), "O");
  edi.InsertAtom(r1, "O3", geom::Vec3(-19.758, -9.617, -18.139), "O");
  edi.InsertAtom(r1, "O4", geom::Vec3(-19.249, -9.772, -20.895), "O");
  edi.InsertAtom(r1, "O5", geom::Vec3(-17.381, -6.775, -19.989), "O");
  edi.InsertAtom(r1, "O6", geom::Vec3(-19.210, -6.695, -22.028), "O");

  // Create mmCIF stream
  writer=MMCifWriter();
  writer.SetStructure(ent, lib, false);
  out = std::stringstream();
  writer.Write("test", out);

  s=out.str();
  // Check that the mmCIF output contains 2 non-polymer entities
  BOOST_CHECK_NE(s.find("loop_\n_entity.id\n_entity.type\n1 non-polymer"),
                 String::npos);
  // Check that atoms are HETATMs since non-poly
  for(auto i: ch.GetAtomList()){
    BOOST_CHECK_NE(s.find("HETATM "+
                          i.GetElement()+" "+
                          i.GetName()+" "+
                          i.GetResidue().GetName()),
                   String::npos);
  }

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();
