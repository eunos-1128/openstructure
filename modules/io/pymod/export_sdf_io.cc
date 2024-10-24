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
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
using namespace boost::python;

#include <ost/export_helper/pair_to_tuple_conv.hh>
#include <ost/io/mol/entity_io_sdf_handler.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/sdf_reader.hh>
#include <ost/io/mol/sdf_writer.hh>
#include <ost/io/sdf_str.hh>
using namespace ost;
using namespace ost::io;
using namespace ost::mol;

String (*sdf_str_a)(const mol::EntityHandle&)=&EntityToSDFString;
String (*sdf_str_b)(const mol::EntityView&)=&EntityToSDFString;

void (*save_sdf_handle)(const mol::EntityHandle& entity, const String& filename)=&SaveSDF;
void (*save_sdf_view)(const mol::EntityView& entity, const String& filename)=&SaveSDF;

void (SDFWriter::*write_handle)(const mol::EntityHandle&)=&SDFWriter::Write;
void (SDFWriter::*write_view)(const mol::EntityView&)=&SDFWriter::Write;

void export_sdf_io()
{
  class_<SDFReader, boost::noncopyable>("SDFReader", init<String, const IOProfile&>())
    .def("Import", &SDFReader::Import)
  ;
  
  class_<SDFWriter, boost::noncopyable>("SDFWriter", init<String>())
    .def("Write", write_handle)
    .def("Write", write_view)
  ;

  // def("LoadSDF", &LoadSDF);
  def("SaveSDF", save_sdf_view);
  def("SaveSDF", save_sdf_handle);

  def("EntityToSDFStr", sdf_str_a);
  def("EntityToSDFStr", sdf_str_b);

  def("SDFStrToEntity", &SDFStringToEntity, (arg("SDF_string"),
                                             arg("profile")=IOProfile()));
}
