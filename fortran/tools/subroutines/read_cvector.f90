! Reads size(A) lines from specificed file, written by write_rvector function

subroutine read_cvector(x, filename, exit_flag)

! Arguments
character(len=*), intent(in) :: filename
complex(dp),intent(inout) :: x(:)
logical, intent(out) :: exit_flag

! Local Variables
real(dp) :: val_re, val_im
integer :: i, u, nx
character(Len=1) :: tab

! exit if file does not exist
INQUIRE(FILE=filename, EXIST=exit_flag)
if( exit_flag .eqv. .False. ) then
    return
end if

! Read Vector
nx = size(x)
u  = 1
open(newunit=u,file=filename, status="old")
do i=1,nx
    read(u,"("//FMT//",a,"//FMT//",a)") val_re, tab, val_im, tab
    x(i) = DCMPLX(val_re, val_im)
enddo
close(u)

end subroutine read_cvector