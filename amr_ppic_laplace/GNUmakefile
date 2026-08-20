# NGA location
NGA_HOME   ?= /home/iliak2/Repositories/Backup/NGA2
FFTW_DIR   = /usr
HYPRE_DIR  = /home/iliak2/Builds/hypre
IRL_DIR    = /home/iliak2/Builds/irl
LAPACK_DIR = /usr
AMREX_DIR  = /home/iliak2/Builds/amrex
HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/openmpi
# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_FFTW  = FALSE
USE_HYPRE = TRUE
USE_LAPACK= FALSE
USE_AMREX = TRUE
USE_HDF5  = TRUE
USE_IRL   = TRUE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga2

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := amrbase amrmpsolvers libraries core
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

# Inform user of Make.packages used
ifdef Ulocs
   $(info Taking user code from: $(Ulocs))
endif
$(info Taking base code from: $(Bdirs))

# Target definition
all: $(executable)
	@echo COMPILATION SUCCESSFUL

# NGA compilation rules
include $(NGA_HOME)/tools/GNUMake/Make.rules
