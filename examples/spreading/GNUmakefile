# NGA location if not yet defined
NGA_HOME ?= ../..

# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_HYPRE = TRUE
USE_LAPACK= TRUE
USE_IRL   = TRUE
USE_ODRPACK= FALSE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# Define external libraries - this can also be in .profile
HYPRE_DIR = /Users/desjardi/Builds/hypre
IRL_DIR   = /Users/desjardi/Builds/IRL
LAPACK_DIR= /Users/desjardi/Builds/lapack
#ODRPACK_DIR=/Users/desjardi/Builds/odrpack

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := core two_phase data solver config grid libraries
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
