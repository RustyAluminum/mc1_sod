subroutine mmc()
    use mShock
    implicit none

    integer i, j                           ! cycle control

    real::dt, dx, t
    real(8):: xrou, xp, xu
    type(pointstate)::allx(1:N)

    character(len=80)::filename = "../output/datamcl.dat"
    open (10, file=filename, status='new')

    ! time&spital size
    dt = 0.01
    dx = 0.1
    t = 0

    !init
    do i = 1, N
        if (i <= 50) then
            xrou = 1.0
            xp = 1.0
            xu = 0
            call init(allx(i), xrou, xp, xu)
        else
            xrou = 0.125
            xp = 0.1
            xu = 0
            call init(allx(i), xrou, xp, xu)
        end if
        write (10, *) t, allx(i)%density, allx(i)%velocity, allx(i)%pressure
    end do

    do while (.true.)
        if (t >= 5) then
            exit
        end if
        call evolutionMC(allx, dx, dt)
        t = t + dt
        do j = 1, N
            write (10, *) t, allx(j)%density, allx(j)%velocity, allx(j)%pressure
        end do
    end do

end subroutine
