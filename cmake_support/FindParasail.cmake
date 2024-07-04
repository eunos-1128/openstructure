# Try to find Parasail
#
# Find the native Parasail include and library. Sets:
# 
#  PARASAIL_INCLUDE_DIR  - Parasail include dir.
#  PARASAIL_LIBRARY      - Parasail library.

find_path(PARASAIL_INCLUDE_DIR parasail.h HINTS PARASAIL_INCLUDE_DIR REQUIRED)
find_library(PARASAIL_LIBRARY NAMES parasail HINTS PARASAIL_LIBRARY REQUIRED)

mark_as_advanced (PARASAIL_LIBRARY PARASAIL_INCLUDE_DIR)
