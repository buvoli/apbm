! Reads size(A) lines from specificed file, written by write_rvector function

subroutine read_rvector(x, filename, exit_flag)

! Arguments
character(len=*), intent(in) :: filename
real(dp),intent(out) :: x(:)
logical, intent(out) :: exit_flag

! Local Variables
real(dp) :: val
integer :: i,u,nx

! exit if file does not exist
INQUIRE(FILE=filename, EXIST=exit_flag)
if( exit_flag .eqv. .False. ) then
    return
end if

! Read Vector
nx = size(x)
open(newunit=u,file=filename, status="old")
do i=1,nx
    read(u,"("//FMT//")") val
    x(i) = val
enddo
close(u)

end subroutine read_rvector