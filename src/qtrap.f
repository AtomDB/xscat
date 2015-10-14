      subroutine qtrap(func, a, b, TotResult, epsilon)
      
      implicit none

      integer MaxStep
      parameter (MaxStep=25)
     
      logical done
      integer iStep
      real Result(MaxStep), error(MaxStep), TotResult
      real a, b, epsilon
      real func
      external func

      done = .false.
      error(1) = 1.0e30
      TotResult = 0.0

      do iStep=1,MaxStep
         if (.not. done) then 
            call trapzd(func, a, b, TotResult, iStep)
            Result(iStep) = TotResult
            if (iStep .gt. 1) then 
               if (Result(iStep-1).ne.0.0) then
                  error(iStep) = abs(Result(iStep)-Result(iStep-1))/
     $                 abs(Result(iStep-1))
c                  if (error(iStep) .gt. error(iStep-1)) then 
c                     write(*,*) 'Error increasing - not good'
c                  end if
               end if
               if (abs(Result(iStep)-Result(iStep-1)) .le.
     $              epsilon*(abs(Result(iStep-1)))) done = .true.
            end if
         end if
      end do

      if (.not. done) then 
         write(*,*) 'No convergence; sequence (step, value, error) is'
         do iStep=1,MaxStep-1
            write(*,*) 'qtrap:',iStep, Result(iStep), error(iStep)
         end do
      end if

      return
      end

      subroutine trapzd(func, a, b, s, n)

      implicit none

      real a, b, s, x
      real*8 del, sum, tnm
      real param(1), r1(1), r2(1), dummy
      integer j, n, it

      real func
      external func

      it = 2**(n-1)
      r1(1) = 0.0
      r2(1) = 0.0

      if (n.eq.1) then
         s = 0.0d0
         param(1) = a
         dummy = func(param, r1, 1)
         param(1) = b
         dummy = func(param, r2, 1)
         s=0.5d0*dble((b-a)*(r1(1) + r2(1)))
      else
         tnm = 1.0d0*it
         del = dble((b-a))/tnm
         x = a + 0.5*real(del)
         sum = 0.0d0
         do j=1,it
c            if (x .gt. b) write(*,*), j, it, x, b
            param(1) = x
            dummy = func(param, r1, 1)
            sum = sum + dble(r1(1))
            x = x+real(del)
         end do
         s = 0.5d0*(s+dble((b-a))*sum/tnm)
      end if
      
      return
      end
