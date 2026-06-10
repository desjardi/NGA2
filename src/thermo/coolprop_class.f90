!> CoolProp-backed material (e.g. IAPWS-95 water) implementing the full class(material) interface.
!> Compiled only when USE_COOLPROP=TRUE. Wraps one low-level CoolProp AbstractState handle per
!> instance via the C bindings in coolprop_cdef, with:
!>   - a cached (input_pair,v1,v2) flash, so the get_p/get_gruneisen pair at one (rho,e) is ONE
!>     CoolProp evaluation instead of two;
!>   - SI<->nondimensional conversion at the boundary (material stores the reference scales, the
!>     solver only ever sees nondimensional quantities). Energy/sound-speed scale on u_ref:
!>     e_ref=u_ref^2, c_ref=u_ref, cv_ref=u_ref^2/T_ref;
!>   - the CLAMP soundness policy: any out-of-domain OR two-phase failure makes get_c return 0,
!>     which the relaxation reads as cL<=0 (unphysical) and bails via ierr. CoolProp signals such
!>     failures through the C errcode (and a _HUGE sentinel), both guarded here -- no NaN reliance.
!>
!> Two-phase liquid policy is CLAMP for now (get_c=0 inside the dome); "project to saturated
!> liquid" is a deliberate TODO, isolated to get_c_from_p_rho / get_gruneisen_from_rho_e.
module coolprop_class
   use precision,        only: WP
   use string,           only: str_medium
   use messager,         only: die
   use material_class,   only: material
   use coolprop_cdef     ! cp_* C bindings + c_long,c_double,c_char,c_null_char
   implicit none
   private

   public :: coolprop

   integer,  parameter :: BUFLEN=512                 !< CoolProp message buffer length
   real(WP), parameter :: HUGE_GUARD=1.0e30_WP       !< reject CoolProp _HUGE sentinel / overflow

   !> Cached AbstractState. Held through a POINTER in coolprop so the const (intent(in)) get_*
   !> methods can still refresh the flash cache (the target is mutable; the association is not).
   type :: cp_state
      integer(c_long) :: handle=-1_c_long
      integer(c_long) :: iDU=0,iPT=0,iDT=0,iDP=0                    !< input-pair keys
      integer(c_long) :: kP=0,kT=0,kD=0,kU=0,kCv=0,kA=0,kH=0,kS=0,kG=0 !< output-parameter keys
      logical         :: cached=.false., last_ok=.false.           !< last-flash cache
      integer(c_long) :: last_pair=-1_c_long
      real(WP)        :: last_v1=0.0_WP, last_v2=0.0_WP
   end type cp_state

   type, extends(material) :: coolprop
      real(WP) :: rho_ref=1.0_WP, u_ref=1.0_WP, T_ref=1.0_WP, p_ref=1.0_WP   !< SI reference scales
      character(len=str_medium) :: backend='BICUBIC&HEOS'   !< default: fastest (tabular over HEOS)
      character(len=str_medium) :: fluid  ='Water'
      type(cp_state), pointer   :: st=>null()
   contains
      procedure :: initialize
      procedure :: get_p_from_rho_e        => cp_p_from_rho_e
      procedure :: get_T_from_p_rho        => cp_T_from_p_rho
      procedure :: get_c_from_p_rho        => cp_c_from_p_rho
      procedure :: get_e_from_p_rho        => cp_e_from_p_rho
      procedure :: get_e_from_p_T          => cp_e_from_p_T
      procedure :: get_p_from_rho_T        => cp_p_from_rho_T
      procedure :: get_rho_from_p_T        => cp_rho_from_p_T
      procedure :: get_cv_from_rho_T       => cp_cv_from_rho_T
      procedure :: get_h_from_p_T          => cp_h_from_p_T
      procedure :: get_hk_from_p_T         => cp_hk_from_p_T
      procedure :: get_s_from_p_T          => cp_s_from_p_T
      procedure :: get_g_from_p_T          => cp_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> cp_gruneisen_from_rho_e
      procedure :: get_T_from_rho_e        => cp_T_from_rho_e
      procedure :: get_c_from_rho_e        => cp_c_from_rho_e
      procedure :: get_cv_from_rho_e       => cp_cv_from_rho_e
      procedure :: get_rhoe_from_p_rho     => cp_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => cp_rhoe_from_p_T
      procedure :: print                   => cp_print
      procedure :: finalize                => cp_finalize
   end type coolprop

contains

   !> Fortran string -> null-terminated C char array
   function cstr(s) result(c)
      character(*), intent(in) :: s
      character(kind=c_char,len=1) :: c(len_trim(s)+1)
      integer :: i
      do i=1,len_trim(s); c(i)=s(i:i); end do
      c(len_trim(s)+1)=c_null_char
   end function cstr

   !> Flash the AbstractState at (input_pair,v1,v2) [SI], reusing the cache if unchanged.
   !> Returns .true. on success; .false. if CoolProp reports an out-of-domain errcode.
   logical function cp_flash(st,pair,v1,v2) result(ok)
      type(cp_state), intent(inout) :: st
      integer(c_long), intent(in)   :: pair
      real(WP),        intent(in)   :: v1,v2
      integer(c_long) :: err
      character(kind=c_char,len=1) :: msg(BUFLEN)
      if (st%cached.and.pair==st%last_pair.and.v1==st%last_v1.and.v2==st%last_v2) then
         ok=st%last_ok; return
      end if
      err=0_c_long
      call cp_update(st%handle,pair,real(v1,c_double),real(v2,c_double),err,msg,int(BUFLEN,c_long))
      ok=(err==0_c_long)
      st%cached=.true.; st%last_pair=pair; st%last_v1=v1; st%last_v2=v2; st%last_ok=ok
   end function cp_flash

   !> Read one keyed output [SI] from the current state; ok=.false. on errcode or _HUGE sentinel.
   real(WP) function cp_out(st,key,ok) result(val)
      type(cp_state), intent(in)  :: st
      integer(c_long), intent(in) :: key
      logical,        intent(out) :: ok
      integer(c_long) :: err
      character(kind=c_char,len=1) :: msg(BUFLEN)
      err=0_c_long
      val=real(cp_keyed_output(st%handle,key,err,msg,int(BUFLEN,c_long)),WP)
      ok=(err==0_c_long.and.abs(val)<HUGE_GUARD)
   end function cp_out

   !> First partial derivative d(Of)/d(Wrt)|Const [SI]; ok=.false. on failure.
   real(WP) function cp_deriv(st,Of,Wrt,Const,ok) result(val)
      type(cp_state), intent(in)  :: st
      integer(c_long), intent(in) :: Of,Wrt,Const
      logical,        intent(out) :: ok
      integer(c_long) :: err
      character(kind=c_char,len=1) :: msg(BUFLEN)
      err=0_c_long
      val=real(cp_first_partial_deriv(st%handle,Of,Wrt,Const,err,msg,int(BUFLEN,c_long)),WP)
      ok=(err==0_c_long.and.abs(val)<HUGE_GUARD)
   end function cp_deriv

   !> Initialize: store reference scales, create the AbstractState handle, cache the integer keys.
   subroutine initialize(this,rho_ref,u_ref,T_ref,p_ref,name,backend,fluid)
      class(coolprop), intent(inout) :: this
      real(WP),        intent(in)    :: rho_ref,u_ref,T_ref,p_ref
      character(*),    intent(in), optional :: name,backend,fluid
      integer(c_long) :: err
      character(kind=c_char,len=1) :: msg(BUFLEN)
      this%rho_ref=rho_ref; this%u_ref=u_ref; this%T_ref=T_ref; this%p_ref=p_ref
      if (present(name))    this%name   =name
      if (present(backend)) this%backend=backend
      if (present(fluid))   this%fluid  =fluid
      this%ns=1
      allocate(this%st)
      err=0_c_long
      this%st%handle=cp_factory(cstr(this%backend),cstr(this%fluid),err,msg,int(BUFLEN,c_long))
      if (err/=0_c_long) call die('[coolprop] AbstractState_factory failed (backend='//trim(this%backend)//', fluid='//trim(this%fluid)//')')
      ! Cache the (constant) input-pair and output-parameter integer keys
      this%st%iDU=cp_input_pair_index(cstr('DmassUmass_INPUTS'))
      this%st%iPT=cp_input_pair_index(cstr('PT_INPUTS'))
      this%st%iDT=cp_input_pair_index(cstr('DmassT_INPUTS'))
      this%st%iDP=cp_input_pair_index(cstr('DmassP_INPUTS'))
      this%st%kP =cp_param_index(cstr('P'))
      this%st%kT =cp_param_index(cstr('T'))
      this%st%kD =cp_param_index(cstr('Dmass'))
      this%st%kU =cp_param_index(cstr('Umass'))
      this%st%kCv=cp_param_index(cstr('Cvmass'))
      this%st%kA =cp_param_index(cstr('speed_of_sound'))
      this%st%kH =cp_param_index(cstr('Hmass'))
      this%st%kS =cp_param_index(cstr('Smass'))
      this%st%kG =cp_param_index(cstr('Gmass'))
   end subroutine initialize

   real(WP) function cp_p_from_rho_e(this,rho,e,y) result(p)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      p=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDU,rho*this%rho_ref,e*this%u_ref**2)) return
      val=cp_out(this%st,this%st%kP,ok); if (ok) p=val/this%p_ref
   end function cp_p_from_rho_e

   !> Temperature from (p,rho): DmassP flash. ILL-CONDITIONED for a near-incompressible liquid
   !> (Halley can fail to converge even at benign states) -- prefer get_T_from_rho_e in any
   !> (rho,e)-driven hot loop. Kept for init/diagnostic callers that genuinely have (p,rho).
   real(WP) function cp_T_from_p_rho(this,p,rho,y) result(T)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      T=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDP,rho*this%rho_ref,p*this%p_ref)) return
      val=cp_out(this%st,this%st%kT,ok); if (ok) T=val/this%T_ref
   end function cp_T_from_p_rho

   !> Temperature from (rho,e): single well-conditioned DmassUmass flash (no DmassP).
   real(WP) function cp_T_from_rho_e(this,rho,e,y) result(T)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      T=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDU,rho*this%rho_ref,e*this%u_ref**2)) return
      val=cp_out(this%st,this%st%kT,ok); if (ok) T=val/this%T_ref
   end function cp_T_from_rho_e

   !> Sound speed from (p,rho): DmassP flash. ILL-CONDITIONED for a stiff liquid -- prefer
   !> get_c_from_rho_e in any (rho,e)-driven hot loop. CLAMP policy: any failure -> 0.
   real(WP) function cp_c_from_p_rho(this,p,rho,y) result(c)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      c=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDP,rho*this%rho_ref,p*this%p_ref)) return
      val=cp_out(this%st,this%st%kA,ok); if (ok) c=val/this%u_ref
   end function cp_c_from_p_rho

   !> Sound speed from (rho,e): single DmassUmass flash. CLAMP: failure (two-phase, etc) -> 0.
   real(WP) function cp_c_from_rho_e(this,rho,e,y) result(c)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      c=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDU,rho*this%rho_ref,e*this%u_ref**2)) return
      val=cp_out(this%st,this%st%kA,ok); if (ok) c=val/this%u_ref
   end function cp_c_from_rho_e

   real(WP) function cp_e_from_p_rho(this,p,rho,y) result(e)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      e=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDP,rho*this%rho_ref,p*this%p_ref)) return
      val=cp_out(this%st,this%st%kU,ok); if (ok) e=val/this%u_ref**2
   end function cp_e_from_p_rho

   real(WP) function cp_e_from_p_T(this,p,T,y) result(e)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      e=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kU,ok); if (ok) e=val/this%u_ref**2
   end function cp_e_from_p_T

   real(WP) function cp_p_from_rho_T(this,rho,T,y) result(p)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      p=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDT,rho*this%rho_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kP,ok); if (ok) p=val/this%p_ref
   end function cp_p_from_rho_T

   real(WP) function cp_rho_from_p_T(this,p,T,y) result(rho)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      rho=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kD,ok); if (ok) rho=val/this%rho_ref
   end function cp_rho_from_p_T

   real(WP) function cp_cv_from_rho_T(this,rho,T,y) result(cv)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      cv=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDT,rho*this%rho_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kCv,ok); if (ok) cv=val*this%T_ref/this%u_ref**2
   end function cp_cv_from_rho_T

   !> cv from (rho,e): single DmassUmass flash (avoids DmassT; uniform (rho,e) routing).
   real(WP) function cp_cv_from_rho_e(this,rho,e,y) result(cv)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      cv=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDU,rho*this%rho_ref,e*this%u_ref**2)) return
      val=cp_out(this%st,this%st%kCv,ok); if (ok) cv=val*this%T_ref/this%u_ref**2
   end function cp_cv_from_rho_e

   real(WP) function cp_h_from_p_T(this,p,T,y) result(h)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      h=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kH,ok); if (ok) h=val/this%u_ref**2
   end function cp_h_from_p_T

   !> Partial specific enthalpies; pure substance (ns=1) -> hk(1)=h.
   subroutine cp_hk_from_p_T(this,p,T,y,hk)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(out) :: hk
      hk=0.0_WP
      hk(1)=this%get_h_from_p_T(p,T,y)
   end subroutine cp_hk_from_p_T

   real(WP) function cp_s_from_p_T(this,p,T,y) result(s)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      s=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kS,ok); if (ok) s=val*this%T_ref/this%u_ref**2
   end function cp_s_from_p_T

   real(WP) function cp_g_from_p_T(this,p,T,y) result(g)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: val; logical :: ok
      g=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      val=cp_out(this%st,this%st%kG,ok); if (ok) g=val/this%u_ref**2
   end function cp_g_from_p_T

   !> Gruneisen Gamma = (1/rho)(dp/de)_rho = (dP/dUmass)|Dmass / rho  (dimensionless).
   real(WP) function cp_gruneisen_from_rho_e(this,rho,e,y) result(G)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: dpdu; logical :: ok
      G=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDU,rho*this%rho_ref,e*this%u_ref**2)) return
      dpdu=cp_deriv(this%st,this%st%kP,this%st%kU,this%st%kD,ok)
      if (ok) G=dpdu/(rho*this%rho_ref)
   end function cp_gruneisen_from_rho_e

   real(WP) function cp_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: u_si; logical :: ok
      rhoe=0.0_WP
      if (.not.cp_flash(this%st,this%st%iDP,rho*this%rho_ref,p*this%p_ref)) return
      u_si=cp_out(this%st,this%st%kU,ok)
      if (ok) rhoe=rho*u_si/this%u_ref**2     ! rho*_nondim * (u_si/e_ref)
   end function cp_rhoe_from_p_rho

   real(WP) function cp_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(coolprop), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho_si,u_si; logical :: ok1,ok2
      rhoe=0.0_WP
      if (.not.cp_flash(this%st,this%st%iPT,p*this%p_ref,T*this%T_ref)) return
      rho_si=cp_out(this%st,this%st%kD,ok1)
      u_si  =cp_out(this%st,this%st%kU,ok2)
      if (ok1.and.ok2) rhoe=rho_si*u_si/(this%rho_ref*this%u_ref**2)
   end function cp_rhoe_from_p_T

   subroutine cp_print(this)
      use,intrinsic :: iso_fortran_env, only: output_unit
      class(coolprop), intent(in) :: this
      write(output_unit,'(a)')        '== coolprop material: '//trim(this%name)
      write(output_unit,'(a)')        '   backend = '//trim(this%backend)//'   fluid = '//trim(this%fluid)
      write(output_unit,'(a,4es12.4)')'   ref(rho,u,T,p) = ',this%rho_ref,this%u_ref,this%T_ref,this%p_ref
   end subroutine cp_print

   subroutine cp_finalize(this)
      class(coolprop), intent(inout) :: this
      integer(c_long) :: err
      character(kind=c_char,len=1) :: msg(BUFLEN)
      if (associated(this%st)) then
         if (this%st%handle>=0_c_long) then
            err=0_c_long
            call cp_free(this%st%handle,err,msg,int(BUFLEN,c_long))
         end if
         deallocate(this%st)
      end if
   end subroutine cp_finalize

end module coolprop_class
