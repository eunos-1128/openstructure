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

#ifndef OST_GDT_HH
#define OST_GDT_HH

#include <ost/base.hh>
#include <ost/geom/geom.hh>

namespace ost { namespace mol { namespace alg {

void GDT(const geom::Vec3List& mdl_pos, const geom::Vec3List& ref_pos,
         int window_size, int max_windows, Real distance_thresh,
         int& n_superposed, geom::Mat4& transform);

}}} // ns

#endif
