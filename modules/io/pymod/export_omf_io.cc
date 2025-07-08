//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2021 by the OpenStructure authors
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

#include <ost/io/mol/omf.hh>
#include <ost/mol/alg/biounit.hh>

using namespace ost;
using namespace ost::io;

namespace{

  template<typename T>
  boost::python::list VecToList(const std::vector<T>& vec){
    boost::python::list l;
    for(typename std::vector<T>::const_iterator it=vec.begin();
        it!=vec.end(); ++it){
      l.append(*it);
    }
    return l;
  }

  PyObject* wrap_to_bytes(OMFPtr omf) {
    String str = omf->ToString();
    return PyBytes_FromStringAndSize(str.c_str(), str.size());
  }
  
  OMFPtr wrap_from_bytes(boost::python::object obj) {
    String str(PyBytes_AsString(obj.ptr()), PyBytes_Size(obj.ptr()));
    return OMF::FromString(str);
  }

  boost::python::list wrap_get_chain_names(OMFPtr omf) {
    return VecToList<String>(omf->GetChainNames());
  }

  boost::python::list wrap_get_chain_types(OMFPtr omf) {

    boost::python::list chain_types;
    std::vector<String> chain_names = omf->GetChainNames();
    for(auto const &cname: chain_names){
      chain_types.append(omf->GetChainData(cname)->chain_type);
    }
    return chain_types;
  }

  OMFPtr wrap_to_assembly(OMFPtr omf,
                          const ost::mol::alg::BUInfo& bu_info) {
    return omf->ToAssembly(bu_info.GetAUChains(),
                           bu_info.GetTransformations(),
                           bu_info.GetBUChains());
  }

  tuple wrap_trace(OMFPtr omf,
                   const String& cname,
                   const String& aname) {
    std::vector<int> rnums;
    geom::Vec3List pos;
    omf->Trace(cname, aname, rnums, pos);
    return boost::python::make_tuple(VecToList(rnums), pos);
  }

  boost::python::list wrap_get_rnums(OMFPtr omf, const String& cname) {
    const std::vector<int>& rnums = omf->GetRNums(cname);
    return VecToList<int>(rnums);
  }

  boost::python::list wrap_get_start_pos(OMFPtr omf, const String& cname) {
    boost::python::list l;
    int i = 0;
    ChainDataPtr chain_data = omf->GetChainData(cname);
    const std::vector<ResidueDefinition>& res_defs = omf->GetResidueDefinitions();
    for(auto res_def_idx: chain_data->res_def_indices) {
      l.append(i);
      i += res_defs[res_def_idx].anames.size();
    }
    return l;
  }
}

void export_omf_io() {

  enum_<OMF::OMFOption>("OMFOption")
    .value("DEFAULT_PEPLIB", OMF::DEFAULT_PEPLIB)
    .value("AVG_BFACTORS", OMF::AVG_BFACTORS)
    .value("ROUND_BFACTORS", OMF::ROUND_BFACTORS)
    .value("SKIP_SS", OMF::SKIP_SS)
    .value("INFER_PEP_BONDS", OMF::INFER_PEP_BONDS)
  ;

  class_<OMF, OMFPtr>("OMF",no_init)
    .def("FromEntity", &OMF::FromEntity, (arg("ent"), arg("max_error")=0.0, arg("options")=0)).staticmethod("FromEntity")
    .def("FromFile", &OMF::FromFile).staticmethod("FromFile")
    .def("FromBytes", &wrap_from_bytes).staticmethod("FromBytes")
    .def("ToFile", &OMF::ToFile)
    .def("ToBytes", &wrap_to_bytes)
    .def("GetMaxError", &OMF::GetMaxError)
    .def("GetAU", &OMF::GetAU)
    .def("GetEntity", &OMF::GetEntity)
    .def("GetAUChain", &OMF::GetAUChain)
    .def("GetEntityChain", &OMF::GetEntityChain)
    .def("GetName", &OMF::GetName)
    .def("GetChainNames", &wrap_get_chain_names)
    .def("GetChainTypes", &wrap_get_chain_types)
    .def("GetPositions", &OMF::GetPositions, return_value_policy<reference_existing_object>(),(arg("cname")))
    .def("GetBFactors", &OMF::GetBFactors, return_value_policy<reference_existing_object>(),(arg("cname")))
    .def("GetOccupancies", &OMF::GetOccupancies, return_value_policy<reference_existing_object>(),(arg("cname")))
    .def("GetResStartPos", &wrap_get_start_pos, (arg("cname")))
    .def("GetAvgBFactors", &OMF::GetAvgBFactors, (arg("cname")))
    .def("GetSequence", &OMF::GetSequence, (arg("cname")))
    .def("GetRNums", &wrap_get_rnums, (arg("cname")))
    .def("ToAssembly", &wrap_to_assembly, (arg("bu_info")))
    .def("Trace", &wrap_trace, (arg("cname"), arg("aname")))
  ;
}
