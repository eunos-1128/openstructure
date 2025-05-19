set(OST_MOL_ALG_PYMOD_SOURCES
  wrap_mol_alg.cc
  export_svd_superpose.cc
  export_clash.cc
  export_trajectory_analysis.cc
  export_structure_analysis.cc
  export_contact_overlap.cc
  export_accessibility.cc
  export_sec_structure.cc
  export_non_standard.cc
  export_molck.cc
  export_membrane.cc
  export_entity_to_density.cc
  export_biounit.cc
  export_gdt.cc
)

set(OST_MOL_ALG_PYMOD_MODULES
  "__init__.py"
  views.py
  superpose.py
  qsscoring.py
  trajectory_analysis.py
  structure_analysis.py
  helix_kinks.py
  hbond.py
  lddt.py
  qsscore.py
  scoring.py
  chain_mapping.py
  stereochemistry.py
  dockq.py
  contact_score.py
  ligand_scoring.py
  ligand_scoring_base.py
  ligand_scoring_scrmsd.py
  ligand_scoring_lddtpli.py
  bb_lddt.py
  scoring_base.py
)

if (NOT ENABLE_STATIC)
  pymod(NAME mol_alg OUTPUT_DIR ost/mol/alg CPP ${OST_MOL_ALG_PYMOD_SOURCES}
        PY ${OST_MOL_ALG_PYMOD_MODULES})
endif()

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure"
                  "stereo_data.json" "STEREO_DATA"
                  "ost_mol_alg_pymod")
install(FILES "stereo_data.json" DESTINATION "share/openstructure/")
