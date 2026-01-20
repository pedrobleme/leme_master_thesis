!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Centralized definition of standard Abaqus integer and real parameters for user subroutines.
!!
!! This module replaces the use of 'vaba_param.inc' and is fully compatible with IMPLICIT NONE, improving code robustness and clarity.
!--------------------------------------------------------------------------------------------------------------------------------------------
      MODULE AbaqusParameters

          IMPLICIT NONE

          ! --- Job control flags ---
          INTEGER, PARAMETER :: J_INT_START_ANALYSIS    = 0
          INTEGER, PARAMETER :: J_INT_START_STEP        = 1
          INTEGER, PARAMETER :: J_INT_SETUP_INCREMENT   = 2
          INTEGER, PARAMETER :: J_INT_START_INCREMENT   = 3
          INTEGER, PARAMETER :: J_INT_END_INCREMENT     = 4
          INTEGER, PARAMETER :: J_INT_END_STEP          = 5
          INTEGER, PARAMETER :: J_INT_END_ANALYSIS      = 6

          ! --- Indices for the integer job information array (i_Array) ---
          INTEGER, PARAMETER :: I_INT_NTOTALNODES    = 1
          INTEGER, PARAMETER :: I_INT_NTOTALELEMENTS = 2
          INTEGER, PARAMETER :: I_INT_KSTEP          = 3
          INTEGER, PARAMETER :: I_INT_KINC           = 4
          INTEGER, PARAMETER :: I_INT_ISTATUS        = 5
          INTEGER, PARAMETER :: I_INT_LWRITERESTART  = 6

          ! --- Indices for the real-valued time information array (r_Array) ---
          INTEGER, PARAMETER :: I_FLT_TOTALTIME = 1
          INTEGER, PARAMETER :: I_FLT_STEPTIME  = 2
          INTEGER, PARAMETER :: I_FLT_DTIME     = 3

      END MODULE AbaqusParameters
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Shared global parameters, flags, and model constants used by Abaqus user subroutines.
!!
!! This module centralizes numerical precision, simulation control flags, physical model parameters, and element-level indexing conventions.
!--------------------------------------------------------------------------------------------------------------------------------------------
      MODULE SharedVariables

          ! Numerical precision
          INTEGER, PARAMETER :: dp = KIND(1.0d0)

          ! Global simulation parameters (set during initialization)
          INTEGER :: n_elems = 0 ! Total number of elements in the mesh (calculated in the subroutine vexternaldb)

          INTEGER, PARAMETER :: max_influence_elems = 500 ! Max elements in nonlocal radius
          INTEGER, PARAMETER :: max_face_neighbors  = 6   ! For C3D8 elements
          INTEGER, PARAMETER :: max_nodes_per_elem  = 8   ! For C3D8 elements

          ! Physical and model parameters
          REAL(KIND=dp), PARAMETER :: nonlocal_length = 0.3_dp ! [mm] Intrinsic length
          REAL(KIND=dp), PARAMETER :: nt = 4000.0 ! Corrosion time in days
          REAL(KIND=dp), PARAMETER :: sigth = 121 ! Material yield stress
          REAL(KIND=dp), PARAMETER :: material_density        = 7.86e-9_dp ! [ton/mm^3] Material density.
          REAL(KIND=dp), PARAMETER :: corroded_volume_limit   = 1000.0_dp    ! Factor for corroded mass criterions [mm^3]
      

          ! Flags for Corrosion Control
          LOGICAL :: failure_occurred_this_increment = .FALSE.
          INTEGER, PARAMETER :: delete_element_flag = 0
          ! States for corrosion control (element_corrosion_status)
          INTEGER, PARAMETER :: flag_corrosion_active = 0
          INTEGER, PARAMETER :: flag_corrosion_inactive = -10

          ! States for strain-based deletion control (strain_deletion_status)
          INTEGER, PARAMETER :: flag_deletion_by_strain_enabled = 0
          INTEGER, PARAMETER :: flag_deletion_by_strain_disabled = -20

          ! States for property update control (property_update_lock)
          INTEGER, PARAMETER :: key_prop_update_unlocked = 0
          INTEGER, PARAMETER :: key_prop_update_locked = -5

          ! Trigger for the global 'key1' variable (General Property Transfer)
          INTEGER, PARAMETER :: flag_general_property_transfer = -8
          INTEGER, PARAMETER :: flag_general_property_transfer_locked = 0

          ! Flag to stop corrosion globally (global_corrosion_stop_flag)
          INTEGER, PARAMETER :: flag_global_stop = -15

          ! Flag for property update
          INTEGER :: prop_update_trigger

          ! Flag to stop global corrosion control
          INTEGER :: global_corrosion_stop_flag

          ! Material Property Parameters
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: elem_properties
          REAL(KIND=dp) :: corroded_mass_total, corrosion_time, corroded_mass_prev

          ! Shared Variables
          INTEGER, DIMENSION(:), ALLOCATABLE :: element_corrosion_status
          INTEGER, DIMENSION(:), ALLOCATABLE :: strain_deletion_status
          INTEGER, DIMENSION(:), ALLOCATABLE :: property_update_lock
          INTEGER, DIMENSION(:), ALLOCATABLE :: n_influence_neighbors
          INTEGER, DIMENSION(:,:), ALLOCATABLE :: influence_map
          REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: influence_distance
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: temp_pitting_write
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: temp_surface_flag_write
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_nonlocal
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_nonlocal_tmp
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_current
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: von_mises_stress
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: non_local_stress
          REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_local

          !-------------------------------------------------------------------
          ! Element property matrix indices (elem_properties)
          !-------------------------------------------------------------------
          INTEGER, PARAMETER :: PROP_ELEM_ID          = 1
          INTEGER, PARAMETER :: PROP_NEIGHBOR_1       = 2
          INTEGER, PARAMETER :: PROP_NEIGHBOR_2       = 3
          INTEGER, PARAMETER :: PROP_NEIGHBOR_3       = 4
          INTEGER, PARAMETER :: PROP_NEIGHBOR_4       = 5
          INTEGER, PARAMETER :: PROP_NEIGHBOR_5       = 6
          INTEGER, PARAMETER :: PROP_NEIGHBOR_6       = 7
          INTEGER, PARAMETER :: PROP_SURFACE_FLAG     = 8
          INTEGER, PARAMETER :: PROP_PITTING_NORM     = 9
          INTEGER, PARAMETER :: PROP_VOLUME           = 10
          INTEGER, PARAMETER :: PROP_COORD_X          = 11
          INTEGER, PARAMETER :: PROP_COORD_Y          = 12
          INTEGER, PARAMETER :: PROP_COORD_Z          = 13
          INTEGER, PARAMETER :: PROP_PITTING_ABS      = 14

      END MODULE SharedVariables
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Abaqus VUMAT entry point.
!!
!! This subroutine acts as a lightweight wrapper for the user-defined material model. Its role is to extract element- and integration-point-
!! specific information from the jblock array and forward all arguments to the core implementation routine (vumatXtrArg).
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vumat (
! Read only -
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

          INCLUDE 'vaba_param.inc'

          DIMENSION jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)

          CHARACTER*80 cmname
          PARAMETER (i_umt_nblock = 1,
     *     i_umt_npt    = 2,
     *     i_umt_layer  = 3,
     *     i_umt_kspt   = 4,
     *     i_umt_noel   = 5 )

      call  vumatXtrArg ( jblock(i_umt_nblock),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(i_umt_noel), jblock(i_umt_npt),
     *     jblock(i_umt_layer), jblock(i_umt_kspt))

      RETURN
      END SUBROUTINE vumat
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Core implementation of the constitutive, corrosion, and damage evolution model.
!!
!! This subroutine contains the full material response algorithm, including plasticity, corrosion-driven degradation, and damage accumulation.
!! It is invoked by the Abaqus VUMAT wrapper after extraction of element- and integration-point-specific data.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vumatXtrArg (
! Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relS  pinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
! Read only extra arguments -
     *     nElement, nMatPoint, nLayer, nSecPoint )

          USE SharedVariables
          INCLUDE 'vaba_param.inc'

!      =========================================================================
!      VUMAT: Variable Declarations
!      =========================================================================
       DIMENSION props(nprops), density(nblock), coordMp(nblock,*),
     1           charLength(nblock), strainInc(nblock,ndir+nshr),
     2           relSpinInc(nblock,nshr), tempOld(nblock),
     3           stretchOld(nblock,ndir+nshr),
     4           defgradOld(nblock,ndir+nshr+nshr),
     5           fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6           stateOld(nblock,nstatev), enerInternOld(nblock),
     7           enerInelasOld(nblock), tempNew(nblock),
     8           stretchNew(nblock,ndir+nshr),
     9           defgradNew(nblock,ndir+nshr+nshr),
     1           fieldNew(nblock,nfieldv),
     2           stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3           enerInternNew(nblock), enerInelasNew(nblock)

!      ! --- Extra Arguments (Indices & Character) ---
       DIMENSION nElement(nblock), nLayer(nblock), nMatPoint(nblock)
       CHARACTER*80 cmname

!      =========================================================================
!      --- Constants & Parameters ---
       REAL(KIND=dp), PARAMETER :: zero   = 0.0_dp
       REAL(KIND=dp), PARAMETER :: third  = 1.0_dp / 3.0_dp
       REAL(KIND=dp), PARAMETER :: half   = 0.50_dp

       ! --- Corrosion Model Parameters ---
       REAL(KIND=dp), PARAMETER :: corrosion_delay_time = 0.0_dp  ! Delay time for corrosion start [days]
       REAL(KIND=dp), PARAMETER :: corrosion_rate       = 0.1417_dp ! Corrosion rate [mm/year]
       ! --- General Damage Model Parameters ---
       REAL(KIND=dp), PARAMETER :: dmax   = 0.999_dp !Maximum damage value before element deletion [dimensionless]
       REAL(KIND=dp), PARAMETER :: beta   = 0.8_dp ! Propagation coefficient

!      =========================================================================
!      --- Local Variables ---

!      --- Elasticity Properties ---
       REAL(KIND=dp) :: E          ! Young's Modulus [Stress]
       REAL(KIND=dp) :: nu         ! Poisson's Ratio [dimensionless]
       REAL(KIND=dp) :: twomu      ! 2 * Shear Modulus (2G) [Stress]
       REAL(KIND=dp) :: alamda     ! Lame's first parameter (λ) [Stress]
       REAL(KIND=dp) :: thremu     ! 3 * Shear Modulus (3G) [Stress]

!      --- Plasticity & Hardening ---
       REAL(KIND=dp) :: yieldOld, yieldNew          ! Yield Stress
       REAL(KIND=dp) :: hard                        ! Hardening modulus (d(sigma_y)/d(epsilon_pl)) [Stress]
       REAL(KIND=dp) :: eq_plastic_strain_inc       ! Increment of PEEQ
       REAL(KIND=dp) :: eq_plastic_strain           ! Accumulated PEEQ
       REAL(KIND=dp) :: peeqOld_k                   ! Old PEEQ (State Var)
       REAL(KIND=dp) :: plasticWorkInc              ! Plastic Work Increment
       REAL(KIND=dp) :: facyld                      ! Yield Factor

!      --- Stress State Variables ---
       REAL(KIND=dp) :: hydrostatic_stress
       REAL(KIND=dp) :: vmises, sigdif
       REAL(KIND=dp) :: trace, stressPower, return_factor
       REAL(KIND=dp) :: s11, s22, s33, s12, s13, s23

!      --- Damage Model (Johnson-Cook & General) ---
       REAL(KIND=dp) :: sigmaJK, epsilonJK          ! JC Stress/Strain terms
       REAL(KIND=dp) :: epsilonf                    ! Failure Strain
       REAL(KIND=dp) :: legacy_damage_old
       REAL(KIND=dp), PARAMETER :: D1 = 0.05_dp  ! Johnson-Cook failure model coefficient D1 [dimensionless]
       REAL(KIND=dp), PARAMETER :: D2 = 3.44_dp  ! Johnson-Cook failure model coefficient D2 [dimensionless]
       REAL(KIND=dp), PARAMETER :: D3 = -2.12_dp ! Johnson-Cook failure model coefficient D3 [dimensionless]
       REAL(KIND=dp), PARAMETER :: D4 = 0.002_dp ! Johnson-Cook failure model coefficient D4 [dimensionless]

!      --- Corrosion & Tribology Model ---
       REAL(KIND=dp) :: corrosion_damage  ! Current corrosion damage state [dimensionless]
       REAL(KIND=dp) :: corrosion_time_days 
       REAL(KIND=dp) :: corr_depth_inc_old
       REAL(KIND=dp) :: damage_tribo_inc
       REAL(KIND=dp) :: masslimit ! General mass-related limit [Mass]
       REAL(KIND=dp), DIMENSION(n_elems) :: corroded_volume

!      --- Geometry & Time Variables ---
       REAL(KIND=dp) :: deltat
       REAL(KIND=dp) :: ku ! Kinetic corrosion parameter [dimensionless]
       REAL(KIND=dp), DIMENSION(n_elems) :: elem_length

!      --- Loop Counters & Status ---
       INTEGER :: k, nvalue, OK


        ! --- DECLARATION OF INDICES FOR STATE VARIABLES (SDVs) ---
  ! SDVs (Solution-Dependent State Variables) are used to store the material state at each integration point.
!-------------------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: SDV_EQ_PLASTIC_STRAIN    = 1  ! Equivalent plastic strain
      INTEGER, PARAMETER :: SDV_TOTAL_DAMAGE         = 2  ! Total accumulated damage
      INTEGER, PARAMETER :: SDV_DELETE_FLAG          = 3  ! Element deletion flag
      INTEGER, PARAMETER :: SDV_PITTING_FACTOR       = 4  ! Pitting factor
      INTEGER, PARAMETER :: SDV_CORROSION_STATUS     = 5  ! Corrosion status (active/inactive)
      INTEGER, PARAMETER :: SDV_NONLOCAL_PROPERTY    = 6  ! Non-local property
      INTEGER, PARAMETER :: SDV_STRESS_S11           = 7  ! Stress component S11
      INTEGER, PARAMETER :: SDV_STRESS_S22           = 8  ! Stress component S22
      INTEGER, PARAMETER :: SDV_STRESS_S33           = 9  ! Stress component S33
      INTEGER, PARAMETER :: SDV_STRESS_S12           = 10 ! Stress component S12
      INTEGER, PARAMETER :: SDV_STRESS_S13           = 11 ! Stress component S13
      INTEGER, PARAMETER :: SDV_STRESS_S23           = 12 ! Stress component S23
      INTEGER, PARAMETER :: SDV_VON_MISES            = 13 ! Von Mises stress
      INTEGER, PARAMETER :: SDV_CORROSION_DEPTH_INC  = 14 ! NEW - V_4_6 Corrosion depth increment
      INTEGER, PARAMETER :: SDV_DSCLd                = 15 ! (Not used, kept for compatibility)
      INTEGER, PARAMETER :: SDV_NONLOCAL_STRESS      = 16 ! Non-local stress
      INTEGER, PARAMETER :: SDV_SURFACE_FLAG         = 17 ! Surface flag
      INTEGER, PARAMETER :: SDV_TOTAL_CORR_DAMAGE    = 18 ! Total corrosion damage
      INTEGER, PARAMETER :: SDV_CURRENT_TOTAL_DAMAGE = 19 ! Current total damage
      INTEGER, PARAMETER :: SDV_JC_DAMAGE_INC        = 20 ! Johnson-Cook damage increment
      INTEGER, PARAMETER :: SDV_CORROSION_DEPTH      = 21 ! Corrosion depth
      INTEGER, PARAMETER :: SDV_CORROSION_START_TIME = 22 ! Corrosion start time
      INTEGER, PARAMETER :: SDV_CORROSION_INIT_FLAG  = 24 ! Corrosion initialization flag
      INTEGER, PARAMETER :: SDV_TRIBO_DAMAGE_ACCUM   = 26 ! Accumulated tribology damage
      INTEGER, PARAMETER :: SDV_CORRODED_VOLUME      = 30 ! Corroded volume

    ! ----------------------------------------------------------------------
    ! Expected parameters in props (elastic properties)
    !   props(1)   - Young's modulus
    !   props(2)   - Poisson's ratio
    ! props(1): Young Modulus, props(2): Poisson's ratio
    ! props(3..) - Hardening curve data points (strain-stress pairs)
  !-----------------------------------------------------------------------
      E = props(1) 
      nu = props(2)
      twomu = E / (1.0_dp + nu)
      alamda = twomu * nu / (1.0_dp - 2.0_dp * nu)
      thremu = 1.5_dp * twomu
      nvalue = (nprops / 2) - 1
      masslimit = material_density * corroded_volume_limit ! Mass limit for corrosion stopping criterion (based on corroded volume limit)

    ! Material behavior logic
      IF (totalTime.EQ.zero) THEN
        ! Initial elastic guess for the first increment - Important for explicit simulations.
        DO k = 1, nblock
            trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            stressNew(k,1) = stressOld(k,1) + twomu * strainInc(k,1) + alamda * trace
            stressNew(k,2) = stressOld(k,2) + twomu * strainInc(k,2) + alamda * trace
            stressNew(k,3) = stressOld(k,3) + twomu * strainInc(k,3) + alamda * trace
            stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
            IF (nshr.GT.1) THEN
                stressNew(k,5) = stressOld(k,5) + twomu * strainInc(k,5)
                stressNew(k,6) = stressOld(k,6) + twomu * strainInc(k,6)
            END IF
        END DO

        corroded_mass_prev = 0.0
       ELSE

        deltat = nt*dt

        ! --- Main Loop over Integration Points ---
        DO k = 1, nblock

            ! Initializes loop variables
            corrosion_damage = 0.0_dp
            damage_tribo_inc = 0.0_dp
            s13 = 0.0_dp
            s23 = 0.0_dp

            peeqOld_k = stateOld(k, SDV_EQ_PLASTIC_STRAIN)
            corr_depth_inc_old = stateOld(k, SDV_CORROSION_DEPTH_INC)
            legacy_damage_old = stateOld(k, SDV_DSCLd)

            damage_local(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)

            ! Corrosion Damage Calculation
            ! The calculation occurs only if the step time is less than 1 (proportional to 1 in Abaqus CAE)
            ! 'flag' is a marker for corrosion calculation control
            IF ((totalTime.LT.1.0_dp) .AND. (element_corrosion_status(nElement(k)) .NE. flag_corrosion_inactive)) THEN
                ! corrosion_time_days: adjusted time considering delay (corrosion_delay_time)
                corrosion_time_days = (totalTime-corrosion_delay_time+dt) * nt
                ! ku: corrosion kinetic parameter (adjustable)
                ku = 1.0
                ! ku = 0.3125 * (corrosion_time_days**(-0.68))  ! Example of experimentally calibrated value

                IF (stateOld(k, SDV_CORROSION_INIT_FLAG) .NE. 10) THEN
                    stateNew(k, SDV_CORROSION_INIT_FLAG) = 10
                    stateNew(k, SDV_CORROSION_START_TIME) = corrosion_time_days
                END IF

                stateNew(k, SDV_CORROSION_DEPTH) = ABS(corrosion_rate *
     1          (corrosion_time_days - stateNew(k, SDV_CORROSION_START_TIME)) / 365.) *       
     2          elem_properties(nElement(k), PROP_SURFACE_FLAG)
     
                IF (stateNew(k, SDV_CORROSION_DEPTH) .LT. stateOld(k, SDV_CORROSION_DEPTH)) THEN
                    stateNew(k, SDV_CORROSION_DEPTH) = stateOld(k, SDV_CORROSION_DEPTH)
                END IF

                elem_length(nElement(k)) = elem_properties(nElement(k), PROP_VOLUME) ** (1./3.)

                corroded_volume(nElement(k)) = ABS(elem_length(nElement(k)) *
     1          elem_length(nElement(k)) * (stateNew(k,SDV_CORROSION_DEPTH))) *
     2          damage_nonlocal(nElement(k))

                IF(corroded_volume(nElement(k)).GE.elem_properties(nElement(k), PROP_VOLUME)) THEN
                    corrosion_damage = dmax
                    stateNew(k,SDV_CORROSION_DEPTH) = 0.0_dp

                ELSE
                    corrosion_damage = (ABS(corroded_volume(nElement(k)))) * ku
                END IF

                damage_local(nElement(k))= corrosion_damage
                corrosion_time = totalTime - corrosion_delay_time ! Accounts for corrosion time (scaling factor)

                !if (abs(corrosion_time*nt) .ge. nt*1.0) then
                ! Interrupts the corrosion process when the corroded mass reaches the defined critical limit (e.g., 100% of total corrosion time).
                ! The block_all routine is called to block corrosion in all elements.

                IF (corroded_mass_total.GE.masslimit) THEN ! original line of code...stop by corroded mass
                    !write(*,*) corroded_mass_total, corrosion_time
                    !write(*,*)'End of Corrosion 100% of volume reached'
                    CALL block_all(ok)
                END IF

                IF (damage_local(nElement(k)).GE.dmax) THEN

                    stateNew(k,SDV_DELETE_FLAG) = delete_element_flag ! Variable that controls the element deletion.
                    damage_local(nElement(k)) = dmax
                    element_corrosion_status(nElement(k)) = flag_corrosion_inactive ! Flag that removes the element from the corrosion calculation.
                    prop_update_trigger = flag_general_property_transfer ! Flag for property update (update only the data of the specific element) when there is element failure
                    CALL process_failed_element_neighbors(nElement(k), beta, OK)
                    !call compute_nonlocal_property(ok) ! Homogeniza os valores locais.
                    failure_occurred_this_increment = .TRUE.
                    stateNew(k,SDV_CORROSION_DEPTH) = 0.0
                END IF
            END IF

            ! Ensures that the damage is not reduced below the previous damage.
            IF (damage_local(nElement(k)).LT.stateOld(k, SDV_TOTAL_DAMAGE)) THEN
                damage_local(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)
            END IF

            ! --- Plasticity Calculation (J2 Mises) ---
            CALL vuhard(yieldOld, hard, peeqOld_k, props(3), nvalue) ! Calculation of yield stress and isotropic hardening.

            ! Elastic trial stress (von Mises yielding criterion)
            trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            s11 = stateOld(k, SDV_STRESS_S11) + twomu*strainInc(k,1) + alamda *trace
            s22 = stateOld(k, SDV_STRESS_S22) + twomu*strainInc(k,2) + alamda*trace
            s33 = stateOld(k, SDV_STRESS_S33) + twomu*strainInc(k,3) + alamda*trace
            s12 = stateOld(k, SDV_STRESS_S12) + twomu * strainInc(k,4)
            IF (nshr.GT.1) THEN
                s13 = stateOld(k, SDV_STRESS_S13) + twomu * strainInc(k,5)
                s23 = stateOld(k, SDV_STRESS_S23) + twomu * strainInc(k,6)
            END IF

            hydrostatic_stress = third * (s11 + s22 + s33)
            s11 = s11 - hydrostatic_stress
            s22 = s22 - hydrostatic_stress
            s33 = s33 - hydrostatic_stress

            IF(nshr .EQ. 1)THEN
                vmises = SQRT(1.5_dp * (s11*s11 + s22*s22 + s33*s33 + 2.0_dp*s12*s12))
            ELSE
                vmises = SQRT(1.5_dp * (s11*s11 + s22*s22 + s33*s33 + 2.0_dp*s12*s12
     1          +2.0_dp*s13*s13 + 2.0_dp*s23*s23))
            END IF

            sigdif = vmises - yieldOld

            IF (sigdif.GT.zero) THEN
                eq_plastic_strain_inc =  sigdif / (thremu + hard) ! Increment of equivalent plastic strain.
            END IF

            ! --- Tribological Damage Calculation (Johnson-Cook) ---
            IF (vmises .GT. 0.0) THEN
                sigmaJK = ((s11 + s22 + s33) / 3.0) / vmises ! Hydrostatic stress / Von Mises
            ELSE
                sigmaJK = 0.0_dp ! If the von Mises stress is zero, there is no tribological damage.
            END IF

            epsilonJK = ABS(trace/3.0) ! epsilon_dot/epsilon_dot_0=1/s^-1
            epsilonf = 0.0_dp
            IF (epsilonJK.NE.0.0) THEN
                epsilonf = (D1 + D2 * EXP(D3 * sigmaJK)) * (1.0 + D4 * LOG(epsilonJK))
                IF (ISNAN(epsilonf))  epsilonf = 0.0
            END IF

            IF ((epsilonf .NE. 0.0) .AND. damage_local(nElement(k)) .LT.dmax) THEN
                damage_tribo_inc = eq_plastic_strain_inc/epsilonf
                damage_local(nElement(k)) = damage_local(nElement(k)) + damage_tribo_inc

                IF (damage_local(nElement(k)) .GE. dmax) THEN
                    stateNew(k,SDV_DELETE_FLAG) = delete_element_flag ! Variable that controls the element deletion.
                    damage_local(nElement(k)) = dmax
                    element_corrosion_status(nElement(k)) = flag_corrosion_inactive
                    prop_update_trigger = flag_general_property_transfer

                    CALL process_failed_element_neighbors(nElement(k), beta, OK)
                    !call compute_nonlocal_property(ok)
                    failure_occurred_this_increment = .TRUE.

                END IF
            END IF

            ! --- Stress Update (Return Mapping + Damage) ---
            yieldNew = yieldOld + hard * eq_plastic_strain_inc
            return_factor= yieldNew / (yieldNew + thremu * eq_plastic_strain_inc)

            stressNew(k,1) = (s11 * return_factor + hydrostatic_stress) * (1 - damage_local(nElement(k)))
            stressNew(k,2) = (s22 * return_factor + hydrostatic_stress) * (1 - damage_local(nElement(k)))
            stressNew(k,3) = (s33 * return_factor + hydrostatic_stress) * (1 - damage_local(nElement(k)))
            stressNew(k,4) = (s12  * return_factor) * (1 - damage_local(nElement(k)))
            IF (nshr .GT. 1) THEN
                stressNew(k,5) = (s13 * return_factor) * (1-damage_local(nElement(k)))
                stressNew(k,6) = (s23 * return_factor) * (1-damage_local(nElement(k)))
            END IF

            eq_plastic_strain = stateOld(k, SDV_EQ_PLASTIC_STRAIN) + eq_plastic_strain_inc
            ! If the total strain exceeds the experimentally obtained value, the element is also removed from the model.
            IF ((eq_plastic_strain .GT. 0.137_dp) .AND. (strain_deletion_status(nElement(k)) .NE.
     1      flag_deletion_by_strain_disabled)) THEN             
              stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
              CALL process_failed_element_neighbors(nElement(k), beta, OK)
              CALL compute_nonlocal_property(ok)
              prop_update_trigger = flag_general_property_transfer
            END IF

            von_mises_stress(nElement(k)) = vmises ! von Mises stress

            ! --- State Variables Update (stateNew) ---
            damage_current(nElement(k)) = damage_local(nElement(k))

            stateNew(k,SDV_EQ_PLASTIC_STRAIN) = stateOld(k,1) + eq_plastic_strain_inc
            stateNew(k,SDV_TOTAL_DAMAGE) = damage_local(nElement(k))
            stateNew(k,SDV_PITTING_FACTOR)=elem_properties(nElement(k), PROP_PITTING_NORM)
            stateNew(k,SDV_CORROSION_STATUS) = element_corrosion_status(nElement(k))
            stateNew(k,SDV_NONLOCAL_PROPERTY) = damage_nonlocal(nElement(k))
            stateNew(k,SDV_STRESS_S11) = (s11 * return_factor + hydrostatic_stress)
            stateNew(k,SDV_STRESS_S22) = (s22 * return_factor + hydrostatic_stress)
            stateNew(k,SDV_STRESS_S33) = (s33 * return_factor + hydrostatic_stress)
            stateNew(k,SDV_STRESS_S12) = (s12 * return_factor)
            stateNew(k,SDV_STRESS_S13) = (s13 * return_factor)
            stateNew(k,SDV_STRESS_S23) = (s23 * return_factor)
            stateNew(k,SDV_VON_MISES) = von_mises_stress(nElement(k))
            stateNew(k,SDV_CORROSION_DEPTH_INC) = corr_depth_inc_old
            stateNew(k,SDV_DSCLd) = legacy_damage_old
            stateNew(k,SDV_NONLOCAL_STRESS) = non_local_stress(nElement(k))
            stateNew(k,SDV_SURFACE_FLAG) = elem_properties(nElement(k), PROP_SURFACE_FLAG)
            stateNew(k,SDV_TOTAL_CORR_DAMAGE) = corrosion_damage
            stateNew(k,SDV_CURRENT_TOTAL_DAMAGE) = damage_current(nElement(k))
            stateNew(k,SDV_JC_DAMAGE_INC) = damage_tribo_inc
            stateNew(k,SDV_CORRODED_VOLUME) = corroded_volume(nElement(k))

            ! --- Energy Update ---

            IF ( nshr .EQ. 1 ) THEN
                stressPower = half * (
     1          ( stateOld(k,SDV_STRESS_S11) + stressNew(k,1) ) * strainInc(k,1) +
     2          ( stateOld(k,SDV_STRESS_S22) + stressNew(k,2) ) * strainInc(k,2) +
     3          ( stateOld(k,SDV_STRESS_S33) + stressNew(k,3) ) * strainInc(k,3) ) +
     4          ( stateOld(k,SDV_STRESS_S12) + stressNew(k,4) ) * strainInc(k,4)
            ELSE
                stressPower = half * (
     1          ( stateOld(k,SDV_STRESS_S11) + stressNew(k,1) ) * strainInc(k,1) +
     2          ( stateOld(k,SDV_STRESS_S22) + stressNew(k,2) ) * strainInc(k,2) +
     3          ( stateOld(k,SDV_STRESS_S33) + stressNew(k,3) ) * strainInc(k,3) ) +
     4          ( stateOld(k,SDV_STRESS_S12) + stressNew(k,4) ) * strainInc(k,4) +
     5          ( stateOld(k,SDV_STRESS_S13) + stressNew(k,5) ) * strainInc(k,5) +
     6          ( stateOld(k,SDV_STRESS_S23) + stressNew(k,6) ) * strainInc(k,6)
            END IF
            enerInternNew(k) = enerInternOld(k) + stressPower / density(k)

            plasticWorkInc = half * ( yieldOld + yieldNew ) * eq_plastic_strain_inc
            enerInelasNew(k) = enerInelasOld(k) + plasticWorkInc / density(k) ! Update of dissipated inelastic energy
        END DO ! End of main loop.
      END IF

      RETURN

      END SUBROUTINE vumatXtrArg
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Calculates the yield stress and isotropic hardening modulus.
!!
!! A interpolation is performed on the yield stress versus equivalent plastic strain curve provided in the material property table.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vuhard(syield, hard, eqplas, table, nvalue)
      USE SharedVariables
      USE AbaqusParameters
      !implicit none
      INCLUDE 'vaba_param.inc'

      ! --- Arguments ---
      INTEGER, INTENT(IN)       :: nvalue
      REAL(KIND=dp), INTENT(IN) :: eqplas, table(2,*)
      REAL(KIND=dp), INTENT(OUT)      :: syield, hard

      ! --- Local Variables ---
      INTEGER :: k1
      REAL(KIND=dp) :: eqpl1, eqpl0, deqpl, syiel0, syiel1, dsyiel
      REAL(KIND=dp), PARAMETER :: zero = 0.0_dp

      ! Initializes output values assuming extrapolation
      syield = table(1, nvalue)
      hard = zero

      ! If the table has more than one point, finds the correct interval.
      IF (nvalue .GT. 1) THEN
          DO k1 = 1, nvalue - 1
              eqpl1 = table(2, k1 + 1)

              IF(eqplas .LT. eqpl1) THEN ! linear interpolation
                  eqpl0 = table(2, k1)
                  syiel0 =table(1, k1)
                  syiel1 = table(1, k1+1)

                  deqpl = eqpl1 - eqpl0
                  dsyiel = syiel1 - syiel0

                  IF (deqpl .GE. 0.0) THEN
                      hard = dsyiel / deqpl
                  ELSE
                      hard = 0.0
                  END IF

                  syield = syiel0  + (eqplas - eqpl0) * hard

                  EXIT
              END IF

          END DO
      END IF
      RETURN
      END SUBROUTINE vuhard
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Entry point for control and data exchange with Abaqus.
!!
!! This subroutine is invoked by Abaqus at predefined analysis events and is used to initialize global variables, read external data, and
! execute control logic.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

      USE SharedVariables
      IMPLICIT NONE

      !include 'vaba_param.inc'
      INTEGER, INTENT(IN) :: lOp, niArray, nrArray, i_Array(niArray)
      REAL(KIND=dp), INTENT(IN) :: r_Array(nrArray)
      INTEGER :: k, p, l, kStep, kInc, Statu, Nel, ok
      INTEGER :: line_count, read_status, unit_num = 100
      CHARACTER(LEN=512) :: line_buffer ! Buffer to read a line

      INTEGER, PARAMETER :: J_INT_START_ANALYSIS    = 0
      INTEGER, PARAMETER :: J_INT_START_STEP        = 1
      INTEGER, PARAMETER :: J_INT_SETUP_INCREMENT   = 2
      INTEGER, PARAMETER :: J_INT_START_INCREMENT   = 3
      INTEGER, PARAMETER :: J_INT_END_INCREMENT     = 4
      INTEGER, PARAMETER :: J_INT_END_STEP          = 5
      INTEGER, PARAMETER :: J_INT_END_ANALYSIS      = 6
    ! ----------------------------------------------------------------------
    ! Indices to access information in the i_Array vector
      INTEGER, PARAMETER :: i_int_nTotalNodes  = 1,  ! Total number of nodes
     1           i_int_nTotalElements = 2,  ! Total number of elements
     2           i_int_kStep          = 3,  ! Current step number
     3           i_int_kInc           = 4,  ! Current increment number
     4           i_int_iStatus        = 5,  ! Analysis status
     5           i_int_lWriteRestart  = 6  ! Flag for restart writing
    ! Possible values for the lOp argument (analysis event)
      INTEGER,PARAMETER ::  j_int_StartAnalysis  = 0,  ! Start of analysis
     1           j_int_StartStep      = 1,  ! Start of step
     2           j_int_SetupIncrement = 2,  ! Setup of increment
     3           j_int_StartIncrement = 3,  ! Start of increment
     4           j_int_EndIncrement   = 4,  ! End of increment
     5           j_int_EndStep        = 5,  ! End of step
     6           j_int_EndAnalysis    = 6  ! End of analysis

    ! Possible values for i_Array(i_int_iStatus)
      INTEGER, PARAMETER ::  j_int_Continue = 0,  ! Continue analysis
     1           j_int_TerminateStep    = 1,  ! Terminate step
     2           j_int_TerminateAnalysis= 2  ! Terminate analysis

    ! Indices to access information in the r_Array vector
      INTEGER, PARAMETER ::  i_flt_TotalTime = 1,  ! Total analysis time
     1           i_flt_StepTime  = 2,  ! Current step time
     2           i_flt_dTime     = 3  ! Time increment

      kStep = i_Array(i_int_kStep)
      kInc  = i_Array(i_int_kInc)
      Statu = i_Array(i_int_iStatus)
      Nel   = i_Array(i_int_nTotalElements)

      SELECT CASE (lOp)

    !---------------------------------------------------------------
    ! Event: Start of Analysis
    !---------------------------------------------------------------
      CASE (j_int_StartAnalysis)


        OPEN(UNIT=UNIT_NUM, FILE='E:\Pedro Bampa\thesis\vumat_big\element_properties.txt',
     1           status='old', action='read', iostat=read_status)


        ! Checks if the file was opened successfully
        IF (read_status .NE. 0) THEN
            ! 1. Prints a highly visible error message in the log.
            PRINT *, '============================================================'
            PRINT *, '>>> FATAL ERROR IN VEXTERNALDB SUBROUTINE <<<'
            PRINT *, '>>> Could not open the properties file Simples.txt.'
            PRINT *, '>>> Check if the file exists in the working directory.'
            PRINT *, '============================================================'

            ! 2. Uses the standard Fortran command to STOP EVERYTHING.
            STOP 'Execution terminated due to file opening error.'
        END IF

        line_count = 0
        DO
            READ(UNIT_NUM, '(A)', IOSTAT=read_status) line_buffer
            ! If iostat > 0 (error) or < 0 (end of file), exit the loop
            IF (read_status /= 0) EXIT
            ! If the line is not blank, count it
            IF (TRIM(line_buffer) /= '') THEN
                line_count = line_count + 1
            END IF
        END DO
        REWIND(UNIT_NUM) ! Rewinds the file to the beginning

        n_elems = line_count
        PRINT *, 'Diagnostic: Properties file loaded. Detected number of elements = ', n_elems

        CALL allocate_arrays(ok)

        DO k = 1, n_elems
            READ(UNIT_NUM, *) (elem_properties(k,l), l=1, PROP_PITTING_ABS)
        END DO
        CLOSE(UNIT_NUM)

        CALL block_surface(ok)
        DO p = 1, n_elems
            IF (elem_properties(p, PROP_SURFACE_FLAG) == 1.0_dp) THEN
                temp_pitting_write(p) = elem_properties(p, PROP_PITTING_NORM)
                temp_surface_flag_write(p) = elem_properties(p, PROP_SURFACE_FLAG)
            ELSE
                elem_properties(p, PROP_PITTING_NORM) = 0.0_dp
                temp_pitting_write(p) = elem_properties(p, PROP_PITTING_NORM)
            END IF
        END DO
        prop_update_trigger = flag_general_property_transfer
        CALL build_influence_map(ok)
        CALL compute_nonlocal_property(ok)
        CALL transfer_properties(ok)

    !---------------------------------------------------------------
    ! Event: end of increment
    !---------------------------------------------------------------
      CASE (j_int_EndIncrement)
        !call transfer_properties(ok)
        !call compute_mass_loss(ok) ! Calculation of mass loss due to corrosion or damage.
        !call compute_nonlocal_stress(ok) ! Calculation of nonlocal stress.

        ! Checks if any failure occurred during the increment
        IF (failure_occurred_this_increment) THEN
            ! 1. Transfers the properties updated by neighbors
            CALL transfer_properties(ok)
            ! 2. Recalculates the nonlocal property for the current state
            CALL compute_nonlocal_property(ok)
            ! 3. Resets the flag for the next increment
            failure_occurred_this_increment = .FALSE.
        END IF

        ! These routines run every increment, regardless of failure
        CALL compute_mass_loss(ok)
        CALL compute_nonlocal_stress(ok)
    !---------------------------------------------------------------
    ! Event: end of analysis
    !---------------------------------------------------------------

      CASE (j_int_EndAnalysis)
          CALL deallocate_arrays(ok)
          CLOSE(100)

      END SELECT
      RETURN
      END SUBROUTINE vexternaldb
  !--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Updates a single property for the neighbors of a failed element.
!!
!! This worker subroutine propagates the effect of a failed element ('nElement') to its direct neighbors. The operation is performed on
!  a specified column of the 'elem_properties' matrix, defined by 'prop_column_index'. The updated values are written to a specified target
! array ('target_write_array').
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE update_neighbor_property(nElement, beta, prop_column_index, target_write_array, ok)
      USE SharedVariables
      IMPLICIT NONE

      ! --- Arguments ---
      INTEGER, INTENT(IN) :: nElement
      REAL(KIND=dp), INTENT(IN) :: beta
      INTEGER, INTENT(IN) :: prop_column_index
      REAL(KIND=dp), INTENT(INOUT) :: target_write_array(n_elems)
      INTEGER, INTENT(OUT) :: ok

      ! --- Local Variables ---
      INTEGER :: j, neighbor_element
      REAL(KIND=dp) :: neighbor_prop_old, prop_from_failed

      ! Gets the property of the failed element to propagate to neighbors
      prop_from_failed = elem_properties(nElement, prop_column_index)

      ! Iterates over the 6 neighbors of the failed element
      DO j = 1, max_face_neighbors
          neighbor_element = elem_properties(nElement, j + 1)

          ! Proceeds only if the neighbor exists (not zero)
          IF (neighbor_element .NE. 0) THEN
              ! If the neighbor is not 'locked', calculates and updates the property
              IF (property_update_lock(neighbor_element) .NE. key_prop_update_locked) THEN
                  neighbor_prop_old = elem_properties(neighbor_element, prop_column_index)

                  ! The logic 'if new < old, use old' is the same as taking the MAX of the two.
                  target_write_array(neighbor_element) = MAX(neighbor_prop_old, prop_from_failed * beta)
              END IF
          END IF
      END DO

      ! Zeros the property in the failed element and sets the locking flag
      target_write_array(nElement) = 0.0_dp
      property_update_lock(nElement) = key_prop_update_locked

      ok = 1
      END SUBROUTINE update_neighbor_property

!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Updates all neighbor associated with a failed element.
!!
!! When an element fails during the VUMAT execution, this routine is invoked to propagate all failure-related effects to its neighboring
!! elements. The procedure delegates each update to the generic worker subroutine 'update_neighbor_property', ensuring consistent treatment
!! of all affected properties (e.g., surface flag and pitting factor).
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE process_failed_element_neighbors(nElement, beta, ok)
      USE SharedVariables
      IMPLICIT NONE

      ! --- Arguments ---
      INTEGER, INTENT(IN) :: nElement
      REAL(KIND=dp), INTENT(IN) :: beta
      INTEGER, INTENT(OUT) :: ok

      ! --- Logic ---
      ! This routine now encapsulates the two actions that always occur together.

      ! 1. Updates the pitting property (column 9), writing to 'temp_pitting_write'.
      CALL update_neighbor_property(nElement, beta, PROP_PITTING_NORM, temp_pitting_write, ok)

      ! 2. Updates the surface flag (column 8), writing to 'temp_surface_flag_write'.
      CALL update_neighbor_property(nElement, beta, PROP_SURFACE_FLAG, temp_surface_flag_write, ok)

      ! Returns OK if both calls are successful (the OK logic can be improved if necessary).
      ok = 1

      END SUBROUTINE process_failed_element_neighbors
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Computes the nonlocal (homogenized) property.
!!
!! Surface elements created due to neighboring element failure receive a homogenized property value based on contributions from surrounding
!! surface elements within a characteristic length.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_property(ok)
      USE SharedVariables
      IMPLICIT NONE

      ! --- Arguments ---
      INTEGER, INTENT(OUT) :: ok

      ! --- Local Variables ---
      INTEGER :: i, p, neighbor_element
      REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
      REAL(KIND=dp) :: dist_sq, nonlocal_length_sq

      nonlocal_length_sq = nonlocal_length**2.0_dp ! Calculates the square once outside the loops

      DO i = 1, n_elems
          ! Proceeds only if element 'i' is surface.
          IF (elem_properties(i, PROP_PITTING_NORM) /= 0.0_dp) THEN

              ! Initializes the sums for element 'i'
              denominator_sum = elem_properties(i, PROP_VOLUME) ! Contribution of the element itself to volume
              numerator_sum = elem_properties(i, PROP_PITTING_NORM)*elem_properties(i, PROP_VOLUME) ! Contribution of the element itself to property
              DO p = 1, n_influence_neighbors(i)
                  neighbor_element = influence_map(i, p)

                  ! Considers only surface neighbors
                  IF (elem_properties(neighbor_element, PROP_PITTING_NORM) /= 0.0_dp) THEN
                      dist_sq = influence_distance(i, p)**2.0_dp
                      weight_p = (1.0_dp - (dist_sq / nonlocal_length_sq))**2.0_dp

                      denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_VOLUME))
                      numerator_sum = numerator_sum + (weight_p * elem_properties(neighbor_element, PROP_PITTING_NORM)
     1                * elem_properties(neighbor_element, PROP_VOLUME))

                  END IF
              END DO

              ! Calculates the final value of the nonlocal property
              IF (denominator_sum > 1.0e-20_dp) THEN
                  damage_nonlocal_tmp(i) = numerator_sum / denominator_sum
              ELSE
                  damage_nonlocal_tmp(i) = 0.0_dp
              END IF

          ELSE
              ! For internal elements, the nonlocal property is zero.
              damage_nonlocal_tmp(i) = 0.0_dp
          END IF
      END DO

      ok = 1

      END SUBROUTINE compute_nonlocal_property
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Transfers updated properties to global storage.
!!
!! This routine copies homogenized and neighbor-updated properties from temporary write arrays into the global property matrices.
!! The transfer is executed only when the global trigger flag is activated.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE transfer_properties(ok)

      USE SharedVariables
      IMPLICIT NONE

      INTEGER :: i
      INTEGER, INTENT(OUT) :: ok

      IF (prop_update_trigger .EQ. flag_general_property_transfer) THEN ! The transfer only occurs if the global 'trigger' is activated.
          ! At the start of the analysis, all elements are assigned prop_update_trigger=-8 to change values in corrosion blocking areas.
          DO i = 1, n_elems
              damage_nonlocal(i) = damage_nonlocal_tmp(i)
              elem_properties(i, PROP_PITTING_NORM) = temp_pitting_write(i)  ! In this part, updates the matwrite values calculated in other routines like the update routine.
              elem_properties(i, PROP_SURFACE_FLAG) = temp_surface_flag_write(i)
          END DO

          ! Deactivates the trigger to avoid unnecessary re-execution.
          ! It will be reactivated when a new element failure occurs
          prop_update_trigger = flag_general_property_transfer_locked
      END IF

      ok = 1

      END SUBROUTINE transfer_properties
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Builds the neighborhood topology required by the nonlocal model by identifying, for each element, all neighboring elements within
!! the prescribed nonlocal interaction radius.
!!
!! This subroutine computes the Euclidean distance between all pairs of elements based on their centroid coordinates and constructs the
!! influence map used by the nonlocal formulation. For each element, neighboring elements located within the intrinsic length scale are
!! stored together with their distances, and the total number of neighbors is recorded. The resulting data structures
!! (influence_map, influence_distance, and n_influence_neighbors) are later used to evaluate nonlocal averages of damage and stress quantities.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE build_influence_map(ok)

      USE SharedVariables
      IMPLICIT NONE
      INTEGER :: t, u, i, j, m
      REAL(KIND=dp) :: aux1
      INTEGER, INTENT(OUT) :: ok

      ! 1. Initializes the influence vectors and matrices.
      DO t = 1, n_elems
          n_influence_neighbors(t) = 0
          DO u = 1, max_influence_elems
              influence_distance(t,u) = 0.0
              influence_map(t,u) = 0
          END DO
      END DO

      ! 2. Calculates the distance between each pair of elements (i, j).
      DO i = 1, n_elems

          m = 1 ! Resets the neighbor n_influence_neighbors for element 'i'

          DO j = 1, n_elems ! loop over the number of elements in the model

              IF (j .NE. i) THEN
                  aux1 = 0.0
                  ! Calculation of the Euclidean distance between all elements
            aux1 = sqrt((elem_properties(j,PROP_COORD_X) - elem_properties(i,PROP_COORD_X))**2.0 +
     1                (elem_properties(j,PROP_COORD_Y) - elem_properties(i, PROP_COORD_Y))**2.0 +
     2                (elem_properties(j,PROP_COORD_Z) - elem_properties(i,PROP_COORD_Z))**2.0)
                ! If element 'j' is within the influence radius, stores it as a neighbor of element 'i'.
                IF (aux1 .LT. nonlocal_length) THEN

                    influence_map(i,m) = j
                    influence_distance(i,m) = aux1
                    n_influence_neighbors(i) = m

                    m = m + 1

                    aux1 = 0.0

                END IF
            END IF
        END DO
      END DO

      ok = 1

      END SUBROUTINE build_influence_map

!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Safely allocates all global shared arrays required by the subroutine.
!! This routine allocates memory for all shared arrays used in the VUMAT subroutine, ensuring that each array is only allocated if it
!! is not already allocated. This prevents runtime errors related to double allocation.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE allocate_arrays(ok)
      USE SharedVariables
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: ok
      INTEGER :: stat ! Variable to check allocation status


      IF (.NOT. ALLOCATED(damage_local))             ALLOCATE(damage_local(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(elem_properties))          ALLOCATE(elem_properties(n_elems, PROP_PITTING_ABS), STAT=stat)
      IF (.NOT. ALLOCATED(damage_nonlocal))          ALLOCATE(damage_nonlocal(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(property_update_lock))     ALLOCATE(property_update_lock(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(influence_map))            ALLOCATE(influence_map(n_elems, max_influence_elems), STAT=stat)
      IF (.NOT. ALLOCATED(influence_distance))       ALLOCATE(influence_distance(n_elems, max_influence_elems), STAT=stat)
      IF (.NOT. ALLOCATED(temp_pitting_write))       ALLOCATE(temp_pitting_write(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(temp_surface_flag_write))  ALLOCATE(temp_surface_flag_write(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(n_influence_neighbors))    ALLOCATE(n_influence_neighbors(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(damage_nonlocal_tmp))      ALLOCATE(damage_nonlocal_tmp(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(damage_current))           ALLOCATE(damage_current(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(element_corrosion_status)) ALLOCATE(element_corrosion_status(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(strain_deletion_status))   ALLOCATE(strain_deletion_status(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(von_mises_stress))         ALLOCATE(von_mises_stress(n_elems), STAT=stat)
      IF (.NOT. ALLOCATED(non_local_stress))         ALLOCATE(non_local_stress(n_elems), STAT=stat)

      ok = 1
      END SUBROUTINE allocate_arrays

!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Releases dynamically allocated global memory at the end of the analysis.
!! This routine deallocates all dynamically allocated arrays in the shared variables module.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE deallocate_arrays(ok)
      USE SharedVariables
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ok

      ! For each matrix, check if it is allocated BEFORE deallocating.
      IF (ALLOCATED(elem_properties))          DEALLOCATE(elem_properties)
      IF (ALLOCATED(damage_nonlocal))          DEALLOCATE(damage_nonlocal)
      IF (ALLOCATED(property_update_lock))     DEALLOCATE(property_update_lock)
      IF (ALLOCATED(influence_map))            DEALLOCATE(influence_map)
      IF (ALLOCATED(influence_distance))       DEALLOCATE(influence_distance)
      IF (ALLOCATED(temp_pitting_write))       DEALLOCATE(temp_pitting_write)
      IF (ALLOCATED(temp_surface_flag_write))  DEALLOCATE(temp_surface_flag_write)
      IF (ALLOCATED(n_influence_neighbors))    DEALLOCATE(n_influence_neighbors)
      IF (ALLOCATED(damage_nonlocal_tmp))      DEALLOCATE(damage_nonlocal_tmp)
      IF (ALLOCATED(damage_current))           DEALLOCATE(damage_current)
      IF (ALLOCATED(element_corrosion_status)) DEALLOCATE(element_corrosion_status)
      IF (ALLOCATED(strain_deletion_status))   DEALLOCATE(strain_deletion_status)
      IF (ALLOCATED(von_mises_stress))         DEALLOCATE(von_mises_stress)
      IF (ALLOCATED(non_local_stress))         DEALLOCATE(non_local_stress)
      IF (ALLOCATED(damage_local))             DEALLOCATE(damage_local)
      ok=1
      END SUBROUTINE deallocate_arrays
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Evaluates the total mass loss due to corrosion at each increment and monitors corrosion progression.
!!
!! This The routine computes the cumulative corroded mass by summing the mass contribution of each element based on its current
!! corrosion damage, material density, and volume.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE compute_mass_loss(ok)

      USE SharedVariables
      IMPLICIT NONE

      REAL(KIND=dp) :: rho, deltam, corroded_mass_inc
      INTEGER :: i, ok

      !---Control Variables---
      rho = 7.86e-09 ! Material density [ton/mm³]
      deltam = 0.0 !(1.06e-010) ! Maximum mass difference between iterations - Stopping criterion.
      corroded_mass_total = 0.0 ! Total corroded mass
      corroded_mass_inc = 0.0 ! Mass variation between iterations

      DO i=1, n_elems
            corroded_mass_total = corroded_mass_total + (damage_current(i) * rho  * elem_properties(i, PROP_VOLUME))
      END DO

      corroded_mass_inc = corroded_mass_total - corroded_mass_prev ! Calculates the mass variation in the increment.
      ! If corrosion has not been stopped globally and the mass change is
      ! significant, prints the status and updates the mass from the previous step.
      IF ((global_corrosion_stop_flag .NE. flag_global_stop).AND.(corroded_mass_inc.GE.deltam)) THEN
          WRITE(*,*) corroded_mass_total, corrosion_time,elem_properties(PROP_ELEM_ID, PROP_VOLUME)
          corroded_mass_prev = corroded_mass_total ! Updates the mass value from the previous step.

      END IF

      ok=1
      END SUBROUTINE compute_mass_loss
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Initializes corrosion activity by selectively blocking elements based on geometric criteria.
!! This subroutine disables corrosion, damage evolution, and strain-based deletion for elements whose centroids lie outside a predefined
!!  geometric region. For such elements, corrosion-related flags are deactivated and property updates are locked.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE block_surface(ok)

      USE SharedVariables
      IMPLICIT NONE

      REAL(KIND=dp) :: coordz1, coordz2
      INTEGER :: k, ok

      ! Geometric control - Removed from the Abaqus model
      coordz1 = 16.8
      coordz2 = -0.2

      DO k=1, n_elems
          ! Checks if the element centroid is outside the corrosion region (<coordx1 or >coordx2), the element will not be corroded
          IF((elem_properties(k, PROP_COORD_Z) .GE. coordz1) .OR. (elem_properties(k, PROP_COORD_Z) .LE. coordz2)) THEN

              ! Corrosion block: Deactivates all relevant flags
              property_update_lock(k) = key_prop_update_locked ! elem_properties will not be updated
              elem_properties(k, PROP_SURFACE_FLAG) = 0.0 ! The surface flag will remain null.
              element_corrosion_status(k) = flag_corrosion_inactive ! The element will not be corroded.
              damage_current(k) = 0.0 ! The element will remain undamaged.
              strain_deletion_status(k) = flag_deletion_by_strain_disabled ! The element will not be deleted, even if the total strain is greater than the limit imposed by the experiment.

          ELSE
              property_update_lock(k) = key_prop_update_unlocked
              element_corrosion_status(k) = flag_corrosion_active
              damage_current(k) = 0.0
              strain_deletion_status(k) = flag_deletion_by_strain_enabled

          END IF
      END DO

      ok=1
      END SUBROUTINE block_surface
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Globally terminates corrosion processes for all elements.
!! When a global stopping criterion is reached, this routine disables corrosion for all elements by locking property updates and setting
!! corrosion status flags to inactive.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE block_all(ok)

      USE SharedVariables
      IMPLICIT NONE
      INTEGER :: k, ok

      DO k=1, n_elems


          property_update_lock(k) = key_prop_update_locked ! flag to not update matwrite
          element_corrosion_status(k) = flag_corrosion_inactive ! flag if equal to -10 does not calculate corrosion


      END DO

      global_corrosion_stop_flag = flag_global_stop ! Activates the global corrosion stop flag.

      ok=1

      END SUBROUTINE block_all
!--------------------------------------------------------------------------------------------------------------------------------------------
!> @brief Computes a nonlocal equivalent stress for surface elements using a weighted neighborhood average.
!! This routine evaluates a nonlocal stress measure for surface elements by averaging the von Mises stresses of neighboring elements within
!!  the nonlocal interaction radius. Internal or inactive elements are excluded from the nonlocal stress calculation.
!--------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_stress(ok)
      USE SharedVariables
      IMPLICIT NONE

      ! --- Arguments ---
      INTEGER, INTENT(OUT) :: ok

      ! --- Local Variables ---
      INTEGER :: i, p, neighbor_element
      REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
      REAL(KIND=dp) :: dist_sq, nonlocal_length_sq

      nonlocal_length_sq = nonlocal_length**2.0_dp

      DO i = 1, n_elems
        ! Considers only unblocked surface elements.
        IF ((elem_properties(i, PROP_PITTING_NORM) /= 0.0_dp) .AND. (element_corrosion_status(i) .NE. flag_corrosion_inactive)) THEN

            ! Initializes the sums for element 'i'
            denominator_sum = elem_properties(i, PROP_VOLUME)
            numerator_sum = von_mises_stress(i) * elem_properties(i, PROP_VOLUME)


            DO p = 1, n_influence_neighbors(i)
                neighbor_element = influence_map(i, p)

                IF (element_corrosion_status(neighbor_element) /= flag_corrosion_inactive) THEN
                    dist_sq = influence_distance(i, p)**2.0_dp
                    weight_p = (1.0_dp - (dist_sq / nonlocal_length_sq))**2.0_dp

                    denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_VOLUME))
                    numerator_sum = numerator_sum + (weight_p * von_mises_stress(neighbor_element) *
     1                elem_properties(neighbor_element, PROP_VOLUME))
                END IF
            END DO

            ! Calculates the final value of the nonlocal stress
            IF (denominator_sum > 1.0e-20_dp) THEN
                non_local_stress(i) = numerator_sum / denominator_sum
            ELSE
                non_local_stress(i) = 0.0_dp
            END IF

            ! Applies the yield stress threshold
            IF (non_local_stress(i) <= sigth) THEN
                non_local_stress(i) = 0.0_dp
            END IF

        ELSE
            ! For internal or inactive elements, the nonlocal stress is zero.
            non_local_stress(i) = 0.0_dp
        END IF
      END DO

      ok = 1
      END SUBROUTINE compute_nonlocal_stress
