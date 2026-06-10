!> Raw ISO_C_BINDING interfaces to the CoolProp C API (CoolPropLib, low-level AbstractState).
!> Compiled only when USE_COOLPROP=TRUE (see src/thermo/Make.package). These are the bare C
!> bindings; coolprop_class wraps them with Fortran-string handling, error checking, a cached
!> (rho,e) flash, and the SI<->nondim conversion. CoolProp must be built with -DEXTERNC so the
!> C API has C linkage (otherwise the symbols are C++-mangled and unlinkable here).
!>
!> All entry points are plain long/double/char* cdecl. char* arguments are C strings: pass a
!> 1-byte character array terminated with C_NULL_CHAR. Every call carries an (errcode,message,
!> buffer_length) trio -- errcode/=0 on failure, with a human-readable reason in message.
module coolprop_cdef
   use,intrinsic :: iso_c_binding, only: c_long,c_double,c_char,c_null_char
   implicit none
   private
   ! Re-export the C kinds and NULL char so coolprop_class needn't use iso_c_binding directly
   public :: c_long,c_double,c_char,c_null_char
   ! C API
   public :: cp_factory,cp_free,cp_update,cp_specify_phase,cp_keyed_output,cp_first_partial_deriv
   public :: cp_input_pair_index,cp_param_index

   interface

      !> AbstractState_factory(backend, fluids, *errcode, message, buffer_length) -> handle
      function cp_factory(backend,fluids,errcode,message,buflen) bind(C,name='AbstractState_factory') result(handle)
         import :: c_long,c_char
         character(kind=c_char,len=1),dimension(*),intent(in)    :: backend,fluids
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
         integer(c_long)                                         :: handle
      end function cp_factory

      !> AbstractState_free(handle, *errcode, message, buffer_length)
      subroutine cp_free(handle,errcode,message,buflen) bind(C,name='AbstractState_free')
         import :: c_long,c_char
         integer(c_long),value                                   :: handle
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
      end subroutine cp_free

      !> AbstractState_update(handle, input_pair, value1, value2, *errcode, message, buffer_length)
      subroutine cp_update(handle,input_pair,v1,v2,errcode,message,buflen) bind(C,name='AbstractState_update')
         import :: c_long,c_double,c_char
         integer(c_long),value                                   :: handle,input_pair
         real(c_double), value                                   :: v1,v2
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
      end subroutine cp_update

      !> AbstractState_specify_phase(handle, phase, *errcode, message, buffer_length)
      subroutine cp_specify_phase(handle,phase,errcode,message,buflen) bind(C,name='AbstractState_specify_phase')
         import :: c_long,c_char
         integer(c_long),value                                   :: handle
         character(kind=c_char,len=1),dimension(*),intent(in)    :: phase
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
      end subroutine cp_specify_phase

      !> AbstractState_keyed_output(handle, param, *errcode, message, buffer_length) -> value
      function cp_keyed_output(handle,param,errcode,message,buflen) bind(C,name='AbstractState_keyed_output') result(val)
         import :: c_long,c_double,c_char
         integer(c_long),value                                   :: handle,param
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
         real(c_double)                                          :: val
      end function cp_keyed_output

      !> AbstractState_first_partial_deriv(handle, Of, Wrt, Constant, *errcode, message, buffer_length) -> value
      function cp_first_partial_deriv(handle,Of,Wrt,Constant,errcode,message,buflen) &
      &        bind(C,name='AbstractState_first_partial_deriv') result(val)
         import :: c_long,c_double,c_char
         integer(c_long),value                                   :: handle,Of,Wrt,Constant
         integer(c_long),                          intent(inout) :: errcode
         character(kind=c_char,len=1),dimension(*),intent(inout) :: message
         integer(c_long),value                                   :: buflen
         real(c_double)                                          :: val
      end function cp_first_partial_deriv

      !> get_input_pair_index(pair_name) -> long key  (e.g. "DmassUmass_INPUTS","PT_INPUTS")
      function cp_input_pair_index(pair) bind(C,name='get_input_pair_index') result(idx)
         import :: c_long,c_char
         character(kind=c_char,len=1),dimension(*),intent(in) :: pair
         integer(c_long)                                      :: idx
      end function cp_input_pair_index

      !> get_param_index(param_name) -> long key  (e.g. "P","T","Dmass","Umass","speed_sound")
      function cp_param_index(param) bind(C,name='get_param_index') result(idx)
         import :: c_long,c_char
         character(kind=c_char,len=1),dimension(*),intent(in) :: param
         integer(c_long)                                      :: idx
      end function cp_param_index

   end interface

end module coolprop_cdef
