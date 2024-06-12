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
#include <ost/bindings/wrap_parasail.hh>
using namespace boost::python;

void export_parasail() {

  def("ParaLocalAlign", &ost::bindings::ParaLocalAlign, (arg("s1"), arg("s2"),
                                                         arg("gap_open_penalty")=11,
                                                         arg("gap_extension_penalty")=1,
                                                         arg("matrix")="blosum62"));

  def("ParaGlobalAlign", &ost::bindings::ParaGlobalAlign, (arg("s1"), arg("s2"),
                                                           arg("gap_open_penalty")=11,
                                                           arg("gap_extension_penalty")=1,
                                                           arg("matrix")="blosum62"));

  def("ParaSemiGlobalAlign", &ost::bindings::ParaSemiGlobalAlign, (arg("s1"), arg("s2"),
                                                                   arg("gap_open_penalty")=11,
                                                                   arg("gap_extension_penalty")=1,
                                                                   arg("matrix")="blosum62"));

  def("ParasailAvailable", &ost::bindings::ParasailAvailable);
}
