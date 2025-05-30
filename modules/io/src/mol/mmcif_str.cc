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
#include <sstream>
#include <ost/io/mol/mmcif_str.hh>
#include <ost/io/mol/mmcif_reader.hh>
#include <ost/io/mmcif_writer.hh>

namespace ost { namespace io {

String EntityToMMCifString(const ost::mol::EntityHandle& ent,
                           const String& data_name, 
                           ost::conop::CompoundLibPtr compound_lib,
                           bool mmcif_conform,
                           const std::vector<MMCifWriterEntity>& entity_info) {
  std::stringstream ss;
  MMCifWriter writer;
  writer.SetStructure(ent, compound_lib, mmcif_conform, entity_info);
  writer.Write(data_name, ss);
  return ss.str();
}

String EntityToMMCifString(const ost::mol::EntityView& ent,
                           const String& data_name,
                           ost::conop::CompoundLibPtr compound_lib,
                           bool mmcif_conform,
                           const std::vector<MMCifWriterEntity>& entity_info) {
  std::stringstream ss;
  MMCifWriter writer;
  writer.SetStructure(ent, compound_lib, mmcif_conform, entity_info);
  writer.Write(data_name, ss);
  return ss.str();
}

String OMFToMMCifString(const ost::io::OMF& omf,
                        const String& data_name,
                        ost::conop::CompoundLibPtr compound_lib,
                        bool mmcif_conform,
                        const std::vector<MMCifWriterEntity>& entity_info) {

  std::stringstream ss;
  MMCifWriter writer;
  writer.SetStructure(omf, compound_lib, mmcif_conform, entity_info);
  writer.Write(data_name, ss);
  return ss.str();
}

std::tuple<mol::EntityHandle, MMCifInfo, ost::seq::SequenceList>
MMCifStringToEntity(const String& mmcif, const IOProfile& profile, bool process) {
  std::stringstream stream(mmcif);
  mol::EntityHandle ent = mol::CreateEntity();
  MMCifReader reader(stream, ent, profile);
  reader.Parse();
  if(profile.processor && process) {
    profile.processor->Process(ent);
  }
  return std::make_tuple(ent, reader.GetInfo(), reader.GetSeqRes());
}

}}
