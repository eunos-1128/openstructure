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
/*
  Author: Tobias Schmidt
 */
#ifndef OST_IO_SDF_READER_HH
#define OST_IO_SDF_READER_HH

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/fstream.hpp>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/io/module_config.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO SDFReader {
public:
  SDFReader(const String& filename);
  SDFReader(const boost::filesystem::path& loc);
  SDFReader(std::istream& instream);

  bool HasNext();

  void Import(mol::EntityHandle& ent);

private:
  void ClearState(const boost::filesystem::path& loc);
  void NextMolecule();

  void ParseAndAddHeader(const String& line, int line_num, mol::EntityHandle& ent,
                         mol::XCSEditor& editor);

  void ParseAndAddAtom(const String& line, int line_num, mol::EntityHandle& ent,
                       bool hetatm, mol::XCSEditor& editor);

  void ParseAndAddBond(const String& line, int line_num, mol::EntityHandle& ent,
                       mol::XCSEditor& editor);

  mol::ChainHandle curr_chain_;
  mol::ResidueHandle curr_residue_;
  int chain_count_;
  int residue_count_;
  int atom_count_;
  int bond_count_;
  int line_num;
  boost::filesystem::ifstream infile_;
  std::istream& instream_;
  boost::iostreams::filtering_stream<boost::iostreams::input>  in_;
};

}}

#endif
