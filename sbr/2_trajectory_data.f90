module trajectory_data
    use kind_module
    implicit none

    type TrajectoryPoint
        !! тип для хранения точки тректории
        real(wp)    :: dland, dcoll, perpn, dalf
        real(wp)    :: vel, tetai
        real(wp)    :: xnpar
        real(wp)    :: rho
        real(wp)    :: vthc, poloidn
        integer     :: izz, iww, jrad

    end type TrajectoryPoint

    integer, parameter :: max_size = 5000

    type Trajectory
        integer size
        !! size of Trajectory
        real(wp) :: tetin  ! tetzap(itr)
        real(wp) :: xmin   ! xmzap(itr)
        real(wp) :: rin    ! rzap(itr)
        real(wp) :: yn3    ! yn3zap(itr)
        real(wp) :: pow    ! powexit
        integer  :: irs    ! irszap(itr)
        integer  :: iw     ! iwzap(itr)
        integer  :: izn    ! iznzap(itr)
        integer  :: mbad
        integer  :: nrefj 
        type(TrajectoryPoint), allocatable ::  points(:)
        !! 
    contains
        procedure :: reset => reset_metod
        procedure :: add_point => add_point_method
        !procedure :: write => write_metod
    end type Trajectory    

    type(Trajectory), pointer :: current_trajectory
contains

end module trajectory_data