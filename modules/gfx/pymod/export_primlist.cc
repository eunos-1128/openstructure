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
using namespace boost::python;

#include <boost/shared_ptr.hpp>

#include <ost/message.hh>
#include <ost/gfx/prim_list.hh>
using namespace ost;
using namespace ost::gfx;

namespace {
  // used numpy support that has been deprecated...
  void add_mesh(PrimList& p, object ova, object ona, object oca, object oia)
  {
    throw Error("AddMesh requires compiled-in numpy support and that has been deprecated");
  }
}

void export_primlist()
{
  class_<PrimList, bases<GfxObj>, boost::shared_ptr<PrimList>, boost::noncopyable>("PrimList", init<const String& >())
    .def("Clear",&PrimList::Clear)
    .def("_add_line",&PrimList::AddLine)
    .def("_add_point",&PrimList::AddPoint)
    .def("_add_sphere",&PrimList::AddSphere)
    .def("_add_cyl",&PrimList::AddCyl)
    .def("_add_text",&PrimList::AddText)
    .def("AddMesh",add_mesh)
    .def("SetColor",&PrimList::SetColor)
    .def("SetDiameter",&PrimList::SetDiameter)
    .def("SetRadius",&PrimList::SetRadius)
  ;

}
