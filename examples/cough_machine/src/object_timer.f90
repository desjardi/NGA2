module object_timer
     use precision,         only: WP
     use string,            only: str_medium,str_long
     implicit none
     private
    
     ! Expose type/constructor/methods
     public :: objtimer

     type :: objtimer
          logical :: amRoot                                                          !< Timer needs to know who's the boss
          character(len=str_medium) :: name='UNNAMED_OBJTIMER'                       !< Name for objtimer
          real(WP) :: vf_wt,vf_wt_total,vf_core_hours                                !< Advance VF WT 
          real(WP) :: lpt_wt,lpt_wt_total,lpt_core_hours                             !< Advance lagrangian particle WT 
          real(WP) :: sgs_wt,sgs_wt_total,sgs_core_hours                             !< Getting sgs dynamic viscosity WT
          real(WP) :: implicit_wt,implicit_wt_total,implicit_core_hours              !< Implicit solver WT
          real(WP) :: pressure_wt,pressure_wt_total,pressure_core_hours              !< Pressure solver WT
          real(WP) :: step_wt,step_wt_total,step_core_hours                          !< WT spent in the domain per time step 
          real(WP), allocatable, dimension(:) :: vf_wt_array                         !< Array to hold WT for vf advance time
          real(WP), allocatable, dimension(:) :: lpt_wt_array                        !< Array to hold WT for lpt advance time
          real(WP), allocatable, dimension(:) :: sgs_wt_array                        !< Array to hold WT for get sgs visc time
          real(WP), allocatable, dimension(:) :: implicit_wt_array                   !< Array to hold WT for implict solver
          real(WP), allocatable, dimension(:) :: pressure_wt_array                   !< Array to hold WT for pressure solver
          real(WP), allocatable, dimension(:) :: step_wt_array                       !< Array to hold WT for domain time steps
     contains
          procedure :: vf_advance_timer                                              !< Method for recording WT advancing VF in the domain
          procedure :: lpt_advance_timer                                             !< Method for recording WT advancing lpt in the domain
          procedure :: sgs_visc_timer                                                !< Method for recording WT getting sgs visc in the domain
          generic   :: implicit_timer=>implicit_timer_tpns,implicit_timer_incomp 
          procedure,private :: implicit_timer_tpns                                   !< Method for recording WT in implicit solver for tpns class
          procedure,private :: implicit_timer_incomp                                 !< Method for recording WT in implicit solver for incomp class
          generic   :: pressure_timer=>pressure_timer_tpns,pressure_timer_incomp
          procedure,private :: pressure_timer_tpns                                   !< Method for recording WT in pressure solver for tpns class
          procedure,private :: pressure_timer_incomp                                 !< Method for recording WT in pressure solver for incomp class
          procedure :: step_timer                                                    !< Method for recording WT during timestep in a domain
     end type objtimer

     !> Declare objtimer constructor
     interface objtimer
          procedure constructor
     end interface objtimer

contains 

     !> Constructor for objtimer
     function constructor(amRoot,name) result(self)
          implicit none
          type(objtimer) :: self
          logical, intent(in) :: amRoot
          character(len=*), optional :: name
          self%amRoot=amRoot
          if (present(name)) self%name=trim(adjustl(name)) 
          ! Timer initiation
          self%vf_wt      =0.0_WP; self%vf_wt_total      =0.0_WP; self%vf_core_hours      =0.0_WP
          self%lpt_wt     =0.0_WP; self%lpt_wt_total     =0.0_WP; self%lpt_core_hours     =0.0_WP
          self%sgs_wt     =0.0_WP; self%sgs_wt_total     =0.0_WP; self%sgs_core_hours     =0.0_WP
          self%implicit_wt=0.0_WP; self%implicit_wt_total=0.0_WP; self%implicit_core_hours=0.0_WP
          self%pressure_wt=0.0_WP; self%pressure_wt_total=0.0_WP; self%pressure_core_hours=0.0_WP
          self%step_wt    =0.0_WP; self%step_wt_total    =0.0_WP; self%step_core_hours    =0.0_WP
          ! Timer array initiation
          allocate(self%vf_wt_array      (1:1)); self%vf_wt_array      =0.0_WP
          allocate(self%lpt_wt_array     (1:1)); self%lpt_wt_array     =0.0_WP
          allocate(self%sgs_wt_array     (1:1)); self%sgs_wt_array     =0.0_WP
          allocate(self%implicit_wt_array(1:1)); self%implicit_wt_array=0.0_WP
          allocate(self%pressure_wt_array(1:1)); self%pressure_wt_array=0.0_WP
          allocate(self%step_wt_array    (1:1)); self%step_wt_array    =0.0_WP
     end function constructor

     subroutine vf_advance_timer(this,cfg,vf)
          use config_class, only: config
          use vfs_class,    only: vfs
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(vfs),      intent(in)    :: vf
          ! Store vf advance wt at current time step
          this%vf_wt=vf%vf_wt
          ! Add current vf advance wt to array
          this%vf_wt_array=[this%vf_wt_array,this%vf_wt]
          ! Sum array to get total time spent on vf advance
          this%vf_wt_total=sum(this%vf_wt_array)
          ! Convert running total time to hours
          this%vf_wt_total=this%vf_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%vf_core_hours=this%vf_wt_total*cfg%nproc
     end subroutine vf_advance_timer

     subroutine lpt_advance_timer(this,cfg,lp)
          use config_class, only: config
          use lpt_class,    only: lpt
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(lpt),      intent(in)    :: lp
          ! Store lpt advance wt at current time step
          this%lpt_wt=lp%lpt_wt
          ! Add current lpt advance wt to array
          this%lpt_wt_array=[this%lpt_wt_array,this%lpt_wt]
          ! Sum array to get total time spent on lpt advance
          this%lpt_wt_total=sum(this%lpt_wt_array)
          ! Convert running total time to hours
          this%lpt_wt_total=this%lpt_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%lpt_core_hours=this%lpt_wt_total*cfg%nproc
     end subroutine lpt_advance_timer

     subroutine sgs_visc_timer(this,cfg,sgs)
          use config_class,   only: config
          use sgsmodel_class, only: sgsmodel
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(sgsmodel), intent(in)    :: sgs
          ! Store sgs get visc wt at current time step
          this%sgs_wt=sgs%sgs_wt
          ! Add current sgs get visc wt to array
          this%sgs_wt_array=[this%sgs_wt_array,this%sgs_wt]
          ! Sum array to get total time spent on sgs get_visc
          this%sgs_wt_total=sum(this%sgs_wt_array)
          ! Convert running total time to hours
          this%sgs_wt_total=this%sgs_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%sgs_core_hours=this%sgs_wt_total*cfg%nproc
     end subroutine sgs_visc_timer

     subroutine implicit_timer_tpns(this,cfg,fs)
          use config_class, only: config
          use tpns_class,   only: tpns
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(tpns),     intent(in)    :: fs
          ! Store solver wt at current time step
          this%implicit_wt=fs%implicit_solver_wt
          ! Add current solver wt to array
          this%implicit_wt_array=[this%implicit_wt_array,this%implicit_wt]
          ! Sum array to get total time spent with solver
          this%implicit_wt_total=sum(this%implicit_wt_array)
          ! Convert running total time to hours
          this%implicit_wt_total=this%implicit_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%implicit_core_hours=this%implicit_wt_total*cfg%nproc
     end subroutine implicit_timer_tpns

     subroutine implicit_timer_incomp(this,cfg,fs)
          use config_class, only: config
          use incomp_class, only: incomp
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(incomp),   intent(in)    :: fs
          ! Store solver wt at current time step
          this%implicit_wt=fs%implicit_solver_wt
          ! Add current solver wt to array
          this%implicit_wt_array=[this%implicit_wt_array,this%implicit_wt]
          ! Sum array to get total time spent with solver
          this%implicit_wt_total=sum(this%implicit_wt_array)
          ! Convert running total time to hours
          this%implicit_wt_total=this%implicit_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%implicit_core_hours=this%implicit_wt_total*cfg%nproc
     end subroutine implicit_timer_incomp

     subroutine pressure_timer_tpns(this,cfg,fs)
          use config_class, only: config
          use tpns_class,   only: tpns
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(tpns),     intent(in)    :: fs
          ! Store solver wt at current time step
          this%pressure_wt=fs%psolv%solver_wt
          ! Add current solver wt to array
          this%pressure_wt_array=[this%pressure_wt_array,this%pressure_wt]
          ! Sum array to get total time spent with solver
          this%pressure_wt_total=sum(this%pressure_wt_array)
          ! Convert running total time to hours
          this%pressure_wt_total=this%pressure_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%pressure_core_hours=this%pressure_wt_total*cfg%nproc
     end subroutine pressure_timer_tpns

     subroutine pressure_timer_incomp(this,cfg,fs)
          use config_class, only: config
          use incomp_class, only: incomp
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          class(incomp),   intent(in)    :: fs
          ! Store solver wt at current time step
          this%pressure_wt=fs%psolv%solver_wt
          ! Add current solver wt to array
          this%pressure_wt_array=[this%pressure_wt_array,this%pressure_wt]
          ! Sum array to get total time spent with solver
          this%pressure_wt_total=sum(this%pressure_wt_array)
          ! Convert running total time to hours
          this%pressure_wt_total=this%pressure_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%pressure_core_hours=this%pressure_wt_total*cfg%nproc
     end subroutine pressure_timer_incomp

     subroutine step_timer(this,cfg)
          use config_class, only: config
          implicit none
          class(objtimer), intent(inout) :: this
          class(config),   intent(in)    :: cfg
          ! Store time step wt
          this%step_wt=cfg%step_wt
          ! Add current time step wt to array
          this%step_wt_array=[this%step_wt_array,this%step_wt]
          ! Sum array to get total time spent in time step
          this%step_wt_total=sum(this%step_wt_array)
          ! Convert running total time to hours
          this%step_wt_total=this%step_wt_total*(1.00_WP/60.00_WP)*(1.00_WP/60.00_WP)
          ! Calculate core-hours
          this%step_core_hours=this%step_wt_total*cfg%nproc
     end subroutine step_timer
          
end module object_timer