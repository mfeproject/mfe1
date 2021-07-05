!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                    !!
!!  Component of MFE1 Version 0.1 -- 15 June 1996     !!
!!  Neil N. Carlson, Dept of Math, Purdue University  !!
!!                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module common_io

  use precision
  use system_io

  implicit none

  integer, public :: input_unit, log_unit, out_unit

  ! Public procedures.
  public  :: read_tagged_data, abort

  private :: read_tagged_rval, read_tagged_rvec, read_tagged_rmtx, &
             read_tagged_ival, read_tagged_ivec

  interface read_tagged_data
     module procedure read_tagged_rval, read_tagged_rvec, read_tagged_rmtx, &
                      read_tagged_ival, read_tagged_ivec
  end interface

  contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  READ_TAGGED_DATA (generic)
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_tagged_rval (data, desc)

      real(wp), intent(out) :: data
      character(*), intent(in), optional :: desc

      integer ios
      character(len=16) tag

      read (input_unit, *, iostat=ios) tag, data
      if (ios /= 0) call read_tagged_error (ios, trim(tag))

      if (present(desc)) &
      write (log_unit, "(t3,a,t32,'=',es12.4)") desc, data

    end subroutine


    subroutine read_tagged_rvec (data, desc)

      real(wp), intent(out) :: data(:)
      character(*), intent(in), optional :: desc

      integer ios
      character(len=16) tag

      read (input_unit, *, iostat=ios) tag, data
      if (ios /= 0) call read_tagged_error (ios, trim(tag))

      if (present(desc)) &
      write (log_unit, "(t3,a,t32,'=',(t33,4es12.4))") desc, data

    end subroutine


    subroutine read_tagged_rmtx (data, desc)

      real(wp), intent(out) :: data(:,:)
      character(*), intent(in), optional :: desc

      integer ios, j
      character(len=16) tag

      read (input_unit, *, iostat=ios) tag, data
      if (ios /= 0) call read_tagged_error (ios, trim(tag))

      if (present(desc)) then
        write (log_unit, "(t3,a,t32,'=',(t33,4es12.4))") desc, data(:,1)
        do j = 2, size(data,2)
          write (log_unit, "(t32,'=',(t33,4es12.4))") data(:,j)
        end do
      end if

    end subroutine


    subroutine read_tagged_ival (data, desc)

      integer, intent(out) :: data
      character(*), intent(in), optional :: desc

      integer ios
      character(len=16) tag

      read (input_unit, *, iostat=ios) tag, data
      if (ios /= 0) call read_tagged_error (ios, trim(tag))

      if (present(desc)) &
      write (log_unit, "(t3,a,t32,'=',i6)") desc, data

    end subroutine


    subroutine read_tagged_ivec (data, desc)

      integer, intent(out) :: data(:)
      character(*), intent(in), optional :: desc

      integer ios
      character(len=16) tag

      read (input_unit, *, iostat=ios) tag, data
      if (ios /= 0) call read_tagged_error (ios, trim(tag))

      if (present(desc)) &
      write (log_unit, "(t3,a,t32,'=',(t34,i5,7(tr1,i5)))") desc, data

    end subroutine


    subroutine read_tagged_error (iostat, tag)

      integer, intent(in) :: iostat
      character(*), intent(in) :: tag

      character(len=32) file

      inquire (unit=input_unit, name=file)

      if (iostat < 0) then

        call abort ( (/ log_unit, stderr /), &
        "End-of-file encountered while reading data from " // trim(file) // ".")

      else

        call abort ( (/ log_unit, stderr /), &
        "Error reading data from " // trim(file) // &
        ".  Part of last line: " // tag // ".")

      end if

    end subroutine read_tagged_error

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !!  ABORT
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine abort (units, mesg)

      integer, intent(in) :: units(:)
      character(*), intent(in) :: mesg

      integer j, unit

      do j = 1, size(units)
        unit = units(j)
        write (unit, '(3a)') "** ", mesg, "  Aborting."
      end do

      stop

    end subroutine abort

end module common_io
