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

    integer, parameter :: max_size = 10000

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

        real(wp) :: tetzap   ! tetzap(itr)
        real(wp) :: xmzap    ! xmzap(itr)
        real(wp) :: rzap     ! rzap(itr)
        real(wp) :: yn3zap   ! yn3zap(itr)
        real(wp) :: powexit  ! powexit
        integer  :: irszap   ! irszap(itr)
        integer  :: iwzap    ! iwzap(itr)
        integer  :: iznzap   ! iznzap(itr)

        type(TrajectoryPoint), allocatable ::  points(:)
        !! 
    contains
        procedure :: init  => init_method
        procedure :: reset => reset_method
        procedure :: add_point => add_point_method
        !procedure :: write => write_metod
    end type Trajectory    

    type(Trajectory), pointer :: current_trajectory
contains
    subroutine init_method(this)
        !! инициализация траетории
        implicit none
        class(Trajectory), intent(inout) :: this
        print *,'инит массива точек:', size(this%points)
        if (allocated(this%points)) deallocate(this%points)
        this%size = 0
        this%nrefj = 0
        this%mbad = 0 
        allocate(this%points(max_size))
    end subroutine

    subroutine reset_method(this, index)
        !! сброс счетчика
        implicit none
        class(Trajectory), intent(inout) :: this
        integer, intent(in) :: index
        this%size = index
    end subroutine

    subroutine add_point_method(this, tpoint)
        !! добавляение новой точнки траектории
        implicit none
        class(Trajectory), intent(inout) :: this
        class(TrajectoryPoint), intent(in) :: tpoint  
        this%size = this%size + 1
        if (this%size > max_size) then
            print *, 'слишком много точек'
            stop
        end if
        this%points(this%size) = tpoint
    end subroutine

end module trajectory_data