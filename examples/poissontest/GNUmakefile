# NGA location if not yet defined
NGA_HOME ?= ../..

# Compilation parameters
PRECISION  = DOUBLE
USE_MPI    = TRUE
USE_FFTW   = TRUE
USE_HYPRE  = TRUE
USE_LAPACK = TRUE
USE_P3DFFT = TRUE
PROFILE    = FALSE
DEBUG      = TRUE
COMP       = gnu
EXEBASE    = nga

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# Define external libraries - this can also be in .profile
#HDF5_DIR =  /opt/hdf5
HYPRE_DIR = /opt/hypre
LAPACK_DIR = /opt/lapack-3.10.1
FFTW_DIR = /opt/fftw-3.3.10
P3DFFT_DIR = /opt/p3dfft-3/build

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := particles core constant_density data solver config grid libraries
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

