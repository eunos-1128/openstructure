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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/gfx/symmetry_node.hh>
using namespace ost;
using namespace ost::gfx;

void export_SymmetryNode()
{
  class_<SymmetryOp>("SymmetryOp", init<const geom::Mat3, const geom::Vec3&>())
    .def(init<const geom::Mat4&>())
  ;
  class_<SymmetryOpList>("SymmetryOpList", init<>())
    .def(vector_indexing_suite<SymmetryOpList>())
  ;
  class_<SymmetryNode, bases<GfxNode>, boost::shared_ptr<SymmetryNode>,
         boost::noncopyable>("SymmetryNode", init<const String&,
                                                   const SymmetryOpList&>())
  ;

}
