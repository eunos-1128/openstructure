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
#ifndef OST_BINDINGS_PARASAIL_H
#define OST_BINDINGS_PARASAIL_H

#include <ost/seq/sequence_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>

namespace ost { namespace seq { namespace alg {

// The following function declarations are intended to mimic the already
// existing function for LocalAligm/GlobalAlign/SemiGlobalAlign
// One example of weirdness is passing the ost::seq::SubstWeightMatrix which
// won't be used in the end. Parasail comes with its own set of matrices
// and the requested one is identified using subst.GetName(). So there is
// plenty of room for optimizations...
ost::seq::AlignmentList ParaLocalAlign(const ost::seq::ConstSequenceHandle& s1,
                                       const ost::seq::ConstSequenceHandle& s2,
                                       ost::seq::alg::SubstWeightMatrixPtr& subst,
                                       int gap_open = -5, int gap_ext = -2);

ost::seq::AlignmentList ParaGlobalAlign(const ost::seq::ConstSequenceHandle& s1,
                                        const ost::seq::ConstSequenceHandle& s2,
                                        ost::seq::alg::SubstWeightMatrixPtr& subst,
                                        int gap_open = -5, int gap_ext = -2);

ost::seq::AlignmentList ParaSemiGlobalAlign(const ost::seq::ConstSequenceHandle& s1,
                                            const ost::seq::ConstSequenceHandle& s2,
                                            ost::seq::alg::SubstWeightMatrixPtr& subst,
                                            int gap_open = -5, int gap_ext = -2);

bool ParasailAvailable();

}}} //ns

#endif
