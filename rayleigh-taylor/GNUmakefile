# NGA location if not yet defined
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    NGA_HOME   ?= /home/iliak2/Repositories/Backup/NGA2
    FFTW_DIR   = /usr
    HYPRE_DIR  = /home/iliak2/Builds/hypre
    IRL_DIR    = /home/iliak2/Desktop/Spline_Workspace/Bezier_Curves/BSpline_Basis_Functions/CppImplementation/IRL/interface-reconstruction-library/install
    LAPACK_DIR = /usr
endif
ifeq ($(UNAME_S),Darwin)
    NGA_HOME   ?= /Users/iliak/Repositories/NGA2
    FFTW_DIR   = /opt/homebrew/opt/fftw
    HYPRE_DIR  = /Users/iliak/Builds/hypre
    IRL_DIR    = /Users/iliak/Builds/irl
    LAPACK_DIR = /opt/homebrew/opt/lapack
endif
# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_FFTW  = TRUE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga

# Directories that contain user-defined code
Udirs   := src ../modules

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# External libraries are defined in .profile/.bashrc/.zshrc, but could be defined here as well

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := core variable_density data transform solver config grid libraries
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
