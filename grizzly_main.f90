!!This is the main program.

program main_grizzly
	use param
	implicit none 
     	include 'mpif.h'
      	integer :: ierr, status(MPI_STATUS_SIZE), runid, groupid, gr_members, n_subgroup, ii, loopmax, tag1, tag2
	integer :: grmin, grmax, i, j, counter

	tag1=1
	tag2=2


      	call MPI_INIT(ierr)
      	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)


if(parameter_space==0) then	!!This is if running for fixed parameter.
	param_zeta= param_zeta_fix
	param_Mmin=param_Mmin_fix
	param_fx= param_fx_fix_1d
	param_al=param_al_fix_1d

	groupid=1
	gr_members=numtasks
	call init_sim()
	if(rank==0)	write(*,*) 'JOB SENT,  PARAMETERS',  param_zeta, param_mmin, param_fx, param_al
	call mpi_grizzly(groupid, gr_members)

else	!!This is if running for parameter space 

	call init_sim()

     	if (rank .eq. 0) then	!creat job arrays
		call read_jobs
	end if

	if(numtasks>gr_members_max) then
		gr_members=gr_members_max
	else
		gr_members=numtasks-1	!!rank 0 will not work	!!group members are the numbers of workers in a sub group
	end if

	n_subgroup = nint(dble(numtasks-1)/dble(gr_members))

	!Distribute the job among other processes

	if(restart_para_est) then
		counter=restart_para_est_id
	else
		counter=0
	end if

	if(manual_terminate) then
		if(nterminate<counter) then
			write(*,*) 'Manual termination is on nterminate<counter'
			stop
		else
			loopmax=int(dble(nterminate-counter)/dble(n_subgroup))+1	!!master node do not contribute to the computation
		end if
	else

		loopmax=int(dble(n_jobs-counter)/dble(n_subgroup))+1	!!master node do not contribute to the computation
	end if


	if(rank==0) then
		write(*,*) 'Number of loop to be run=', loopmax
		write(*,*) 'Number of sub-group, Members of each group, total cores', n_subgroup, gr_members, numtasks

	end if




	DO ii=1, loopmax
		if(rank ==0) then
			write(*,*) 'Running loop # of #', ii, loopmax
			DO i=1, n_subgroup
				counter=counter+1
				param_zeta=job_arr(counter,1)
				param_mmin=job_arr(counter,2)
				param_mmin=10.d0**param_mmin
				param_fx=job_arr(counter,3)
				param_al=job_arr(counter,4)

				write(*,*) 'JOB SENT, SUB-GROUP, COUNTER, PARAMETERS', i, counter, param_zeta, param_mmin, param_fx, param_al

				groupid=i
				grmin=(i-1)*gr_members+1
				grmax=i*gr_members
				if(grmax>numtasks-1) grmax=numtasks-1
				DO j=grmin, grmax
				call MPI_SEND(counter, 1, MPI_INTEGER, j, tag1, MPI_COMM_WORLD, ierr)
				call MPI_SEND(groupid, 1, MPI_INTEGER, j, tag1, MPI_COMM_WORLD, ierr)
				call MPI_SEND(gr_members, 1, MPI_INTEGER, j, tag1, MPI_COMM_WORLD, ierr)
				call MPI_SEND(param_zeta,1,MPI_DOUBLE_PRECISION,j,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(param_mmin,1,MPI_DOUBLE_PRECISION,j,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(param_fx,1,MPI_DOUBLE_PRECISION,j,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(param_al,1,MPI_DOUBLE_PRECISION,j,tag1,MPI_COMM_WORLD,ierr)
				END DO

			END DO

		
			DO i=1, numtasks-1
				call MPI_RECV(counter, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
			END DO 

			IF(counter>n_jobs) then
				write(*,*) 'EXISTING LOOP'
				!exit
			END IF


		else
			call MPI_RECV(counter, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(groupid, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(gr_members, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(param_zeta, 1, MPI_DOUBLE_PRECISION, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(param_mmin, 1, MPI_DOUBLE_PRECISION, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(param_fx, 1, MPI_DOUBLE_PRECISION, 0, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(param_al, 1, MPI_DOUBLE_PRECISION, 0, tag1, MPI_COMM_WORLD, status, ierr)

			IF(counter <= n_jobs) THEN
				call mpi_grizzly(groupid, gr_members)
			END IF

			call MPI_SEND(counter, 1, MPI_INTEGER, 0, tag2, MPI_COMM_WORLD, ierr)

		end if
	END DO

end if
	if(rank==0)	write(*,*) 'PROGRAM ENDS'
      call MPI_FINALIZE(ierr)

end program main_grizzly




!!This one create the job array and store in common array job_arr
!!Paramters are M, z, delta, fx, alpha..time will be set in rt code.
SUBROUTINE read_jobs
	use param
	IMPLICIT NONE
	integer:: i,j,k,l,m,con



	call read_parameter(zeta_points, nn_zeta, zeta_arr)
	call read_parameter(mmin_points, nn_mmin, mmin_arr)
	call read_parameter(fx_points, nn_fx, fx_arr)
	call read_parameter(al_points, nn_al, al_arr)

	con=0
	DO i=1,nn_zeta
	DO j=1,nn_mmin
	DO k=1,nn_fx
	DO l=1,nn_al
		con=con+1

		job_arr(con,1)=zeta_arr(i)
		job_arr(con,2)=mmin_arr(j)
		job_arr(con,3)=fx_arr(k)
		job_arr(con,4)=al_arr(l)

	END DO
	END DO
	END DO
	END DO

END SUBROUTINE read_jobs

!!This reads parameter files
!
subroutine read_parameter(arr2d, n, arr)
	use param
	implicit none
	real(8), dimension(2), intent(in)::arr2d
	integer, intent(in)::n
	real(8), dimension(n), intent(out)::arr
	integer::fstat,i
	real(8)::del

	if(n>1) then
		del=(arr2d(2)-arr2d(1))/dble(n-1)
	else
		del=0.d0
	end if


	DO i=1, n
	arr(i)=arr2d(1)+dble(i-1)*del
	END DO

	write(*,*) 'parameter space'

	DO i=1,n
	    write(*,'(f20.1)') arr(i)
	END DO

end subroutine read_parameter






