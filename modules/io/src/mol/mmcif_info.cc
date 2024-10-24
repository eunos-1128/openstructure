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

#include <ost/io/io_exception.hh>
#include <ost/io/mol/mmcif_info.hh>
#include <ost/log.hh>

namespace ost { namespace io {

void MMCifInfo::AddMMCifPDBChainTr(String cif, String pdb, bool fault_tolerant)
{
  std::map<String, String>::iterator tr_it = cif_2_pdb_chain_id_.find(cif);
  if (tr_it != cif_2_pdb_chain_id_.end()) {
    std::stringstream msg;
    msg << "mmCIF chain id '" <<  cif << "' is already mapped to '" <<
           tr_it->second << "'. Cannot map it to '" <<  pdb << "'." ;
    if(fault_tolerant) {
      LOG_WARNING(msg.str());
      return;
    } else {
      throw IOException(msg.str());
    }
  } else {
    cif_2_pdb_chain_id_.insert(std::pair<String, String>(cif, pdb));
  }
}

String MMCifInfo::GetMMCifPDBChainTr(String cif) const
{
  std::map<String, String>::const_iterator tr_it =
    cif_2_pdb_chain_id_.find(cif);
  if (tr_it == cif_2_pdb_chain_id_.end()) { return ""; }
  return tr_it->second;
}

void MMCifInfo::AddPDBMMCifChainTr(String pdb, String cif, bool fault_tolerant)
{
  std::map<String, String>::iterator tr_it = pdb_2_cif_chain_id_.find(pdb);
  if (tr_it != pdb_2_cif_chain_id_.end()) {
    std::stringstream msg;
    msg << "PDB chain id '" <<  pdb << "' is already mapped to '" <<
           tr_it->second << "'. Cannot map it to '" <<  cif << "'." ;
    if(fault_tolerant) {
      LOG_WARNING(msg.str());
      return;
    } else {
      throw IOException(msg.str());
    }
  } else {
    pdb_2_cif_chain_id_.insert(std::pair<String, String>(pdb, cif));
  }
}

String MMCifInfo::GetPDBMMCifChainTr(String pdb) const
{
  std::map<String, String>::const_iterator tr_it =
    pdb_2_cif_chain_id_.find(pdb);
  if (tr_it == pdb_2_cif_chain_id_.end()) { return ""; }
  return tr_it->second;
}

void MMCifInfo::AddMMCifEntityIdTr(String cif, String ent_id, bool fault_tolerant)
{
  std::map<String, String>::iterator tr_it = cif_2_entity_id_.find(cif);
  if (tr_it != cif_2_entity_id_.end()) {
    std::stringstream msg;
    msg << "mmCIF chain id '" <<  cif << "' is already mapped to " <<
           "entity id '" << tr_it->second << "'. Cannot map it to '" <<
           ent_id << "'." ;
    if(fault_tolerant) {
      LOG_WARNING(msg.str());
      return;
    } else {
      throw IOException(msg.str());
    }
  } else {
    cif_2_entity_id_.insert(std::pair<String, String>(cif, ent_id));
  }
}

String MMCifInfo::GetMMCifEntityIdTr(String cif) const
{
  std::map<String, String>::const_iterator tr_it =
    cif_2_entity_id_.find(cif);
  if (tr_it == cif_2_entity_id_.end()) { return ""; }
  return tr_it->second;
}

void MMCifInfo::AddAuthorsToCitation(StringRef id, std::vector<String> list,
                                     bool fault_tolerant)
{
  // find citation
  std::vector<MMCifInfoCitation>::iterator cit_it;
  for (cit_it = citations_.begin(); cit_it != citations_.end(); ++cit_it) {
    String cit_id = cit_it->GetID(); // to ensure lifetime of StringRef-pointers
    if (id == StringRef(cit_id.c_str(), cit_id.length())) {
      cit_it->SetAuthorList(list);
      return;
    }
  }

  if(fault_tolerant) {
    LOG_WARNING("No citation for identifier '" + id.str() + "' found. "
                "Couldn't set author list.");
  } else {
    throw IOException("No citation for identifier '" + id.str() + "' found.");
  }
}

void MMCifInfo::AddBioUnit(MMCifInfoBioUnit bu)
{
  std::vector<MMCifInfoBioUnit>::iterator bu_it;
  for (bu_it = biounits_.begin(); bu_it != biounits_.end(); ++bu_it) {
    if (bu_it->GetID() == bu.GetID()) {
      break;
    }
  }
  if (bu_it != biounits_.end()) {
    bu_it->Merge(bu);
    return;
  }
  biounits_.push_back(bu);
}

void MMCifInfoBioUnit::AddChain(String chain)
{
  chains_.push_back(chain);

  if (tr_chains_.size()) {
    tr_chains_.back().second = chains_.size();
  }
  else {
    std::pair<int, int> tr_interval = std::pair<int, int>(0, 1);
    tr_chains_.push_back(tr_interval);
  }
}

void MMCifInfoBioUnit::SetChainList(std::vector<String> chains)
{
  chains_ = chains;

  if (tr_chains_.size()) {
    tr_chains_.clear();
  }
  std::pair<int, int> tr_interval = std::pair<int, int>(0, chains_.size());
  tr_chains_.push_back(tr_interval);
}

void MMCifInfoBioUnit::AddOperations(std::vector<MMCifInfoTransOpPtr> operations)
{
  operations_.push_back(operations);

  if (tr_operations_.size()) {
    tr_operations_.back().second = operations_.size();
  }
  else {
    std::pair<int, int> tr_interval = std::pair<int, int>(0, 1);
    tr_operations_.push_back(tr_interval);
  }
}

void MMCifInfoBioUnit::Merge(MMCifInfoBioUnit& from)
{
  // merge chains
  int old_size = chains_.size();
  chains_.insert(chains_.end(), from.chains_.begin(), from.chains_.end());
  std::pair<int, int> tr_interval = std::pair<int, int>(old_size,
                                                        chains_.size());
  tr_chains_.push_back(tr_interval);
  // merge operations
  old_size = operations_.size();
  operations_.insert(operations_.end(),
                     from.operations_.begin(),
                     from.operations_.end());
  tr_interval.first = old_size;
  tr_interval.second = operations_.size();
  tr_operations_.push_back(tr_interval);
}

MMCifInfoStructRefSeqPtr 
MMCifInfoStructRef::AddAlignedSeq(const String& aid, const String& chain_name, 
                                  int seq_begin, int seq_end, int db_begin, 
                                  int db_end, bool fault_tolerant)
{
  std::map<String, MMCifInfoStructRefSeqPtr>::const_iterator i=seqs_.find(aid);
  if (i!=seqs_.end()) {
    std::stringstream msg;
    msg << "Duplicate align_id for struct_ref '" << id_ << "'";
    if (fault_tolerant) {
      msg << ". Record will be overwritten.";
      LOG_WARNING(msg.str());
    }
    else {
      throw IOException(msg.str());
    }

  }
  MMCifInfoStructRefSeqPtr p(new MMCifInfoStructRefSeq(aid, chain_name,
                                                       seq_begin, seq_end, 
                                                       db_begin, db_end));
  seqs_[aid]=p;
  return p;
}


MMCifInfoStructRefSeqPtr 
MMCifInfoStructRef::GetAlignedSeq(const String& aid) const
{

  std::map<String, MMCifInfoStructRefSeqPtr>::const_iterator i=seqs_.find(aid);
  return i==seqs_.end() ? MMCifInfoStructRefSeqPtr() : i->second;
}

MMCifInfoStructRefSeqDifPtr 
MMCifInfoStructRefSeq::AddDif(int seq_rnum, const String& db_rnum, const String& details)
{
  MMCifInfoStructRefSeqDifPtr d(new MMCifInfoStructRefSeqDif(seq_rnum, db_rnum,
                                                             details));
  difs_.push_back(d);
  return d;
}

void MMCifInfo::AddEntityBranchLink(String chain_name,
                                    int rnum1, int rnum2,
                                    const String& aname1,
                                    const String& aname2,
                                    unsigned char bond_order) {
  // [] operator creates new value if no such element exists for key
  entity_branches_[chain_name].push_back(MMCifInfoEntityBranchLink(rnum1, rnum2,
                                                                   aname1, aname2,
                                                                   bond_order));
}

const std::vector<MMCifInfoEntityBranchLink> MMCifInfo::GetEntityBranchByChain(
                                                 const String& chain_name) const
{
  // search the requested chain
  MMCifInfoEntityBranchLinkMap::const_iterator blm_it =
                                               entity_branches_.find(chain_name);
  if (blm_it != entity_branches_.end()) {
    return blm_it->second;
  }
  return std::vector<MMCifInfoEntityBranchLink>();
}

const std::vector<String> MMCifInfo::GetEntityBranchChainNames() const
{
  std::vector<String> chain_names;
  MMCifInfoEntityBranchLinkMap::const_iterator blm_it;
  for (blm_it = entity_branches_.begin();
       blm_it != entity_branches_.end(); ++blm_it) {
    chain_names.push_back(blm_it->first);
  }
  return chain_names;
}

const MMCifEntityDesc& MMCifInfo::GetEntityDesc(const String& entity_id) const {
  MMCifEntityDescMap::const_iterator it = entity_desc_.find(entity_id);
  if(it == entity_desc_.end()) {
    throw IOException("No EntityDesc for entity id \""+entity_id+"\"");
  }
  return it->second;
}

void MMCifInfo::SetEntityDesc(const String& entity_id,
                              const MMCifEntityDesc& entity_desc)
{
  if(entity_desc_.find(entity_id) != entity_desc_.end()) {
    throw IOException("EntityDesc for entity_id \""+entity_id+"\" already set");
  }
  entity_desc_[entity_id] = entity_desc;
}

std::vector<String> MMCifInfo::GetEntityIds() const {
  std::vector<String> vec;
  for(auto it: entity_desc_) {
    vec.push_back(it.first);
  }
  return vec;
}

std::vector<String> MMCifInfo::GetEntityIdsOfType(const String& entity_type) const {
  std::vector<String> vec;
  for(auto it: entity_desc_) {
    if(it.second.entity_type == entity_type) {
      vec.push_back(it.first);
    }
  }
  return vec;
}

std::ostream& operator<<(std::ostream& os, const MMCifInfoEntityBranchLink& eb)
{
  os << "<MMCifInfoEntityBranchLink rnum1:" << eb.rnum1 <<
     " rnum2: " << eb.rnum2 << " aname1: " << eb.aname1 <<
     " aname2: " << eb.aname2 << "bond_order: " << eb.bond_order << '>';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const std::vector<MMCifInfoEntityBranchLink>& eb_list)
{
  os << "<MMCifInfoEntityBranchLinkList";
  std::vector<MMCifInfoEntityBranchLink>::const_iterator bl_it;
  for (bl_it = eb_list.begin(); bl_it != eb_list.end(); ++bl_it) {
    os << *bl_it;
  }
  os << ">";
  return os;
}

}} //ns
