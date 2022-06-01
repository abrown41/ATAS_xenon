!! gaussian :: produce gaussian shaped, two colour pulse for read in to RMT
!! compilation :: gfortran -o gaussian.x gaussian.f90
!! usage :: ./gaussian.x <delay between pulses in femtoseconds>
!!
module precisn
implicit none
public
    INTEGER, PARAMETER    :: wp = SELECTED_REAL_KIND(12)  ! `double' precision

end module precisn

program gaussian
use precisn, only: wp
implicit none
integer, parameter :: nsteps = 350000
real(wp) :: field(nsteps)
real(wp) :: time(nsteps), shift_time(nsteps)
real(wp) :: IR(nsteps), XUV(nsteps)
real(wp) :: dt = 0.01, midtime, FWHM
real(wp) :: IR_intensity, intensity , E0, frequency
real(wp) :: delay_in_fs, delay_in_au 
character(6) :: delayinp, intensityinp
integer :: ii

delay_in_fs = 1.0_wp

call get_command_argument(1, delayinp)
read(delayinp, *) delay_in_fs
delay_in_au = delay_in_fs * 27.212_wp

call get_command_argument(2, intensityinp)
read(intensityinp, *) IR_intensity

field = 0.0_wp
do ii = 1, nsteps
    time(ii) = (dt/2.0_wp) +  (ii-1) * dt
end do

! IR
midtime = 900.0_wp
FWHM = 155.0_wp
frequency = 0.06798_wp
E0 = 0.05336_wp * sqrt(IR_intensity)
shift_time = time - midtime - delay_in_au
IR = E0 * exp( - (shift_time * shift_time) / (FWHM*FWHM)) * cos(frequency*shift_time)


! XUV
FWHM = 7.0_wp
intensity = 0.20_wp
frequency = 2.05_wp
E0 = 0.05336_wp * sqrt(intensity)
shift_time = time - midtime
XUV = E0 * exp( - (shift_time * shift_time) / (FWHM*FWHM)) * cos(frequency*shift_time)

open (unit=11,file='EField.inp',status='new')
do ii=1, nsteps
    write(11, *) time(ii), 0.0_wp, 0.0_wp, IR(ii) + XUV(ii)
end do
close(11)
end program gaussian
