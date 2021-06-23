module mShock
    implicit none

    type pointstate
        real(8):: density, pressure, velocity
        real(8) :: c, rouu
    end type

    integer, parameter:: N = 100
    real, parameter::gamma = 1.4                    ! gas adiabatic index

contains
    subroutine init(x, xrou, xp, xu)
        type(pointstate)::x
        real(8):: xrou, xp, xu
        x%density = xrou
        x%pressure = xp
        x%velocity = xu

        x%c = sqrt(gamma*x%pressure/x%density)
        x%rouu = xrou*xu
    end subroutine

    subroutine evolutionMC(x, delta_x, delta_t)
        real ::delta_x, delta_t
        type(pointstate)::x(:)
        integer i

        real(8):: temp1_density(N), temp1_velocity(N), temp1_pressure(N)
        real(8):: temp2_density, temp2_velocity, temp2_pressure
        real(8)::tempc1, tempc2
        real(8):: a
        a = delta_t/delta_x

        do i = 2, N - 1
            !step pre
            temp1_density(i) = x(i)%density - a*(2*x(i + 1)%rouu - 2*x(i)%rouu)
            temp1_velocity(i) = x(i)%velocity - a*((x(i + 1)%velocity**2 + x(i + 1)%pressure/x(i + 1)%density) &
                                                   - (x(i)%velocity**2 + x(i)%pressure/x(i)%density))
            temp1_pressure(i) = x(i)%pressure - a* &
                                ( &
                                (x(i + 1)%density*x(i + 1)%velocity*x(i + 1)%c*x(i + 1)%c + x(i + 1)%velocity*x(i + 1)%pressure) &
                                - (x(i + 1)%density*x(i)%velocity*x(i)%c*x(i)%c + x(i)%velocity*x(i)%pressure) &
                                )
        end do

        temp1_density(1) = x(1)%density; temp1_velocity(1) = x(1)%velocity; temp1_pressure(1) = x(1)%pressure; 
        temp1_density(N) = x(N)%density; temp1_velocity(N) = x(N)%velocity; temp1_pressure(N) = x(N)%pressure; 
        do i = 2, N - 1
            !step correction
            !i
            tempc1 = sqrt(gamma*x(i)%pressure/x(i)%density)
            !i-1
            tempc2 = sqrt(gamma*x(i - 1)%pressure/x(i - 1)%density)

            temp2_density = x(i)%density - a*(2*temp1_density(i)*temp1_velocity(i) &
                                              - 2*temp1_density(i - 1)*temp1_velocity(i - 1))
            temp2_velocity = x(i)%velocity - a*((temp1_velocity(i)**2 + temp1_pressure(i)/temp1_density(i)) &
                                                - (temp1_velocity(i - 1)**2 + temp1_pressure(i - 1)/temp1_density(i - 1)))
            temp2_pressure = x(i)%pressure - a* &
                             ( &
                             (temp1_density(i)*temp1_velocity(i)*tempc1*tempc1 + temp1_velocity(i)*temp1_pressure(i)) &
                        - (temp1_density(i - 1)*temp1_velocity(i - 1)*tempc2*tempc2 + temp1_velocity(i - 1)*temp1_pressure(i - 1)) &
                             )

            !step result
            x(i)%density = 0.5*(x(i)%density + temp2_density)
            x(i)%velocity = 0.5*(x(i)%velocity + temp2_velocity)
            x(i)%pressure = 0.5*(x(i)%pressure + temp2_pressure)

            if (x(i)%pressure <= 0) then
                x(i)%pressure = 1.e-50
            end if
        end do

        x(1) = x(2)
        x(N) = x(N - 1)
    end subroutine

end module mShock
