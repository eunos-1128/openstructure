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


#include <boost/python.hpp>

#include <ost/mol/alg/gdt.hh>

using namespace boost::python;

namespace{

boost::python::tuple WrapGDT(const geom::Vec3List& mdl_pos,
                             const geom::Vec3List& ref_pos,
                             int window_size, int max_windows,
                             Real distance_thresh) {
  int n_superposed;
  geom::Mat4 transform;
  ost::mol::alg::GDT(mdl_pos, ref_pos, window_size, max_windows,
                     distance_thresh, n_superposed, transform);
  return boost::python::make_tuple(n_superposed, transform);
}

} // ns

void export_GDT() {

    def("GDT", &WrapGDT, (arg("mdl_pos"), 
                          arg("ref_pos"),
                          arg("window_size"),
                          arg("max_windows"),
                          arg("distance_thresh")));
}
