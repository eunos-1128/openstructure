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

using namespace ost;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

conop::CompoundLibPtr SetDefaultCompoundLib() {
  // return NULL if not successful, else return newly set default lib
  // REQ: OST_ROOT to be set
  char * ost_root = getenv("OST_ROOT");
  if (!ost_root) return conop::CompoundLibPtr();
  SetPrefixPath(ost_root);
  String lib_path = GetSharedDataPath() + "/compounds.chemlib";
  conop::CompoundLibPtr compound_lib = conop::CompoundLib::Load(lib_path);
  if (compound_lib) {
    conop::Conopology::Instance().SetDefaultLib(compound_lib);
  }
  return compound_lib;
}

BOOST_AUTO_TEST_CASE(mmcif_writer_entity1)
{
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
  mol::AtomHandle a=edi.InsertAtom(r, "CA", geom::Vec3(32.0, -128.0, -2.5),
                                   "C", 1.0, 128.0);
  edi.SetChainType(ch, mol::CHAINTYPE_UNKNOWN);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, SetDefaultCompoundLib(), false);
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
  BOOST_TEST_MESSAGE("  Running mmcif_writer_poly_vs_non_poly tests...");
  /*
    Go for small polymers that are not polymer... the story of 2 amino acids (to
    be handled like RCSB as 2 separated chains) plus how the same thing works
    for nucleic acids and branched molecules (sugars).
   */

  // Polypeptide: 2aa in a chain are 2 separated entities for RCSB (check 1E8K)
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor edi=ent.EditXCS();
  mol::ChainHandle ch=edi.InsertChain("A");
  edi.SetChainType(ch, mol::CHAINTYPE_POLY_PEPTIDE_L);
  // ALA
  mol::ResidueHandle r1=edi.AppendResidue(ch, "ALA");  // is_hetatm=false
  edi.InsertAtom(r1, "N",  geom::Vec3(44.987, 17.389, 12.362), "N", 0.81,
                 25.57);
  edi.InsertAtom(r1, "CA", geom::Vec3(45.936, 16.434, 12.890), "C", 0.81,
                 28.21);
  edi.InsertAtom(r1, "C",  geom::Vec3(47.196, 17.227, 13.152), "C", 0.81,
                 33.78);
  edi.InsertAtom(r1, "O",  geom::Vec3(47.506, 18.153, 12.401), "O", 0.81,
                 23.02);
  edi.InsertAtom(r1, "CB", geom::Vec3(46.244, 15.293, 11.961), "C", 0.81,
                 29.85);
  // PRO
  mol::ResidueHandle r2=edi.AppendResidue(ch, "PRO");  // is_hetatm=false
  edi.InsertAtom(r2, "N",   geom::Vec3(47.953, 16.910, 14.229), "N", 0.81,
                 32.19 );
  edi.InsertAtom(r2, "CA",  geom::Vec3(47.673, 15.829, 15.187), "C", 0.81,
                 33.89 );
  edi.InsertAtom(r2, "C",   geom::Vec3(46.564, 16.052, 16.233), "C", 0.81,
                 37.97 );
  edi.InsertAtom(r2, "O",   geom::Vec3(46.059, 17.169, 16.417), "O", 0.81,
                 34.43 );
  edi.InsertAtom(r2, "CB",  geom::Vec3(49.054, 15.755, 15.880), "C", 0.81,
                 36.60 );
  edi.InsertAtom(r2, "CG",  geom::Vec3(49.357, 17.213, 16.030), "C", 0.81,
                 34.77 );
  edi.InsertAtom(r2, "CD",  geom::Vec3(49.098, 17.714, 14.637), "C", 0.81,
                 34.62 );
  edi.InsertAtom(r2, "OXT", geom::Vec3(46.144, 15.129, 16.950), "O", 0.81,
                 34.03 );
  conop::HeuristicProcessor heu_proc;
  // Make sure that the two residues r1, r2 are actually connected
  heu_proc.Process(ent);
  BOOST_CHECK_EQUAL(mol::InSequence(r1, r2), true);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, SetDefaultCompoundLib(), false);
  std::stringstream out;
  writer.Write("test", out);

  // Check that the mmCIF output contains 2 non-polymer entities
  String s=out.str();
  BOOST_CHECK_NE(
        s.find("loop_\n_entity.id\n_entity.type\n1 non-polymer\n2 non-polymer"),
        String::npos);
  BOOST_CHECK_NE(
             s.find("loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2"),
             String::npos);
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

  // write mmCIF, check atoms are marked HETATM

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();
