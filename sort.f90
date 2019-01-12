	SUBROUTINE sort(n, arr)
	IMPLICIT NONE
	integer::n
	REAL(4), DIMENSION(n), INTENT(INOUT) :: arr
	INTEGER, PARAMETER :: NN=15, NSTACK=50
	REAL(4) :: a
	INTEGER :: k,i,j,jstack,l,r
	INTEGER, DIMENSION(NSTACK) :: istack

	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				a=arr(j)
				do i=j-1,l,-1
					if (arr(i) <= a) exit
					arr(i+1)=arr(i)
				end do
				arr(i+1)=a
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap_r(arr(k),arr(l+1))
			call masked_swap_rs(arr(l),arr(r),arr(l)>arr(r))
			call masked_swap_rs(arr(l+1),arr(r),arr(l+1)>arr(r))
			call masked_swap_rs(arr(l),arr(l+1),arr(l)>arr(l+1))
			i=l+1
			j=r
			a=arr(l+1)
			do
				do
					i=i+1
					if (arr(i) >= a) exit
				end do
				do
					j=j-1
					if (arr(j) <= a) exit
				end do
				if (j < i) exit
				call swap_r(arr(i),arr(j))
			end do
			arr(l+1)=arr(j)
			arr(j)=a
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	END SUBROUTINE sort

	SUBROUTINE swap_r(a,b)
	!use nrtype
	implicit none
	REAL(4), INTENT(INOUT) :: a,b
	REAL(4) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap_r

	SUBROUTINE masked_swap_rs(a,b,mask)
	!use nrtype
	implicit none
	REAL(4), INTENT(INOUT) :: a,b
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	LOGICAL(LGT), INTENT(IN) :: mask
	REAL(4) :: swp
	if (mask) then
		swp=a
		a=b
		b=swp
	end if
	END SUBROUTINE masked_swap_rs

!BL
	SUBROUTINE nrerror(string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	write (*,*) 'nrerror: ',string
	STOP 'program terminated by nrerror'
	END SUBROUTINE nrerror
