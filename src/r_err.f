
      subroutine r_err(dumpcore,errstring)

      integer dumpcore
      external abort
      character errstring*(*)

      if (dumpcore.eq.1) then
         write(*,*) 'An error occured, and the program is now'
         write(*,*) 'going to bail out.  The error was:'
         write(*,*) errstring
c         call abort()
      end if
      
      if (dumpcore.eq.0) then
         write(*,*) 'Message: ',errstring
      end if

      return
      end 
