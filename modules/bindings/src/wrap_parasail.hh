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

namespace ost { namespace bindings {

ost::seq::AlignmentHandle ParaLocalAlign(const ost::seq::SequenceHandle& s1,
                                         const ost::seq::SequenceHandle& s2,
                                         int gap_open_penalty = 11,
                                         int gap_extension_penalty = 1,
                                         const String& matrix="blosum62");

ost::seq::AlignmentHandle ParaGlobalAlign(const ost::seq::SequenceHandle& s1,
                                          const ost::seq::SequenceHandle& s2,
                                          int gap_open_penalty = 11,
                                          int gap_extension_penalty = 1,
                                          const String& matrix="blosum62");

ost::seq::AlignmentHandle ParaSemiGlobalAlign(const ost::seq::SequenceHandle& s1,
                                              const ost::seq::SequenceHandle& s2,
                                              int gap_open_penalty = 11,
                                              int gap_extension_penalty = 1,
                                              const String& matrix="blosum62");

bool ParasailAvailable();

}} //ns

#endif
