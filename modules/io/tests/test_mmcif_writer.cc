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

#include <ost/platform.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/mol.hh>
#include <ost/io/mol/mmcif_writer.hh>

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
                                   "C");
  a.SetOccupancy(1.0);
  a.SetBFactor(128.0);
  edi.SetChainType(ch, mol::CHAINTYPE_UNKNOWN);

  // Create mmCIF stream
  MMCifWriter writer;
  writer.SetStructure(ent, SetDefaultCompoundLib(), false);
  std::stringstream out;
  writer.Write("test", out);

  // Check if entity starts with 1, either by reading mmCIF or "grep"
  String s=out.str();
  BOOST_CHECK_EQUAL(
                    s.substr(0, 53),
                    "data_test\nloop_\n_entity.id\n_entity.type\n1 non-polymer"
                    );

  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();
