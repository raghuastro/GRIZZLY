!!This determines the group among the called processors 
!! to do a job under head processor 'leader'
!!
subroutine det_leader(gpid, npgp, leader, gpext)
	use param
	implicit none
	integer, intent(in) :: gpid, npgp
	integer, intent(out) :: leader, gpext

	if(parameter_space==0) then
		leader = (gpid-1)*npgp 
		gpext = leader + npgp-1
	else
		leader = (gpid-1)*npgp +1
		gpext = leader + npgp-1
	end if

end subroutine det_leader


!!This is the main subroutine. Distribute and do all the jobs.
!!
subroutine mpi_grizzly(gpid, npgp)
	use param
	implicit none 
     	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
      	integer :: ierr, status(MPI_STATUS_SIZE), runid, exit_id
	real(8) :: z, zcomm, mf

	call init_grizzly(gpid, npgp)


	runid=1
	exit_id=0

	DO redshift_loop_index=check_min,num_checkpoints
		if(redshift_loop_index==num_checkpoints) exit_id=1

		if(selected_reion==1) then	!!This is if tracking history is not included
			z=z_selectpoint(runid)
			if(track_age==1) then
				call part1(gpid, npgp)
				call part2(gpid, npgp)
				if(z==comm_redshift) call part3(gpid, npgp)
				if(z==comm_redshift) call part4(gpid, npgp)
				if(z==comm_redshift) call part5(gpid, npgp)
				if(z==comm_redshift) call part6(gpid, npgp)
				call part7(gpid, npgp)
				if(z==comm_redshift) runid=runid+1
				if(runid>num_selectpoints) exit_id=1
				call part8(gpid, npgp, exit_id)
			else
				if(z==z_checkpoint(redshift_loop_index)) then
				call part1(gpid, npgp)
				call part2(gpid, npgp)
				call part3(gpid, npgp)
				call part4(gpid, npgp)
				call part5(gpid, npgp)
				call part6(gpid, npgp)
				call part7(gpid, npgp)
				runid=runid+1
				if(runid>num_selectpoints) exit_id=1
				call part8(gpid, npgp, exit_id)
				end if
			end if
			if(runid>num_selectpoints) exit
			if(exit_id==1) exit
		else
				call part1(gpid, npgp)
				call part2(gpid, npgp)
				call part3(gpid, npgp)
				call part4(gpid, npgp)
				call part5(gpid, npgp)
				call part6(gpid, npgp)
				call part7(gpid, npgp)
				call part8(gpid, npgp, exit_id)
				if(exit_id==1) exit
		end if

	END DO
end subroutine mpi_grizzly

!!initialization of the simulation
!
subroutine init_sim()	
	use param
	implicit none
	include 'mpif.h'
	integer ::  ierr,  status(MPI_STATUS_SIZE)
	logical, parameter :: restart = .false.
	integer, parameter :: restart_id=33


      	if (rank .eq. 0) then
		call sim_specification
		call meminfo
		call readcube
		call read_redshift
		if(selected_reion==1) call read_select
		call cal_age_timestep

		if(genxhimap==1 .or. genalphamap==1) then
			haloagenew=0
			halomassnew=0.0
		end if

		if(selected_reion==1) call MPI_BCAST(num_selectpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		if(selected_reion==1) call MPI_BCAST(z_selectpoint,num_selectpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(z_checkpoint,num_checkpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(age_checkpoint,num_checkpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(spdim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	else
		if(selected_reion==1) call MPI_BCAST(num_selectpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		if(selected_reion==1) call MPI_BCAST(z_selectpoint,num_selectpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(num_checkpoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(z_checkpoint,num_checkpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(age_checkpoint,num_checkpoints,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(spdim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(cube,max_cell_read*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	end if	!rank

	if(restart) then
		check_min=restart_id
		!call read_restartarrays(check_min, z_checkpoint(check_min))
	else
		check_min=1
	end if	!restart
end subroutine init_sim


!!initialization of the simulation
!
subroutine init_grizzly(gpid, npgp)	!!gpid is group id, npgp is number of workers in a group
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
	integer ::i, tag1, ierr,  status(MPI_STATUS_SIZE)

	call det_leader(gpid, npgp, leader, gpext)
	if(rank==leader) write(*,*) 'groupid, leader, workers', gpid, leader, npgp 

	tag1=1

      	if (rank .eq. leader) then
		if(useanalytic==0 .or. gentkmap==1) call read_ifront(rank)
		call det_filename
	end if

     	if (rank .eq. leader) then
		DO i=leader+1,gpext
			call MPI_SEND(filename_comm, 180, MPI_CHARACTER, i,tag1, MPI_COMM_WORLD,ierr)
			if(useanalytic==0 .or. gentkmap==1) then
				call MPI_SEND(n_ifront,1,MPI_INTEGER,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(ifront_data,n_ifront*nparam1d,MPI_REAL,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(age1d_max,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(age1d_min,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(mass1d_max,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(mass1d_min,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(delta1d_max,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(delta1d_min,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(z1d_max,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
				call MPI_SEND(z1d_min,1,MPI_DOUBLE_PRECISION,i,tag1,MPI_COMM_WORLD,ierr)
			end if
		END DO
	elseif(rank>leader .and. rank<=gpext) then
		call MPI_RECV(filename_comm, 180, MPI_CHARACTER, leader, tag1, MPI_COMM_WORLD, status, ierr)
		if(useanalytic==0 .or. gentkmap==1 ) then
			call MPI_RECV(n_ifront,1,MPI_INTEGER,leader, tag1, MPI_COMM_WORLD, status, ierr)
			allocate(ifront_data(n_ifront, nparam1d))
			call MPI_RECV(ifront_data,n_ifront*nparam1d,MPI_REAL,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(age1d_max,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(age1d_min,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(mass1d_max,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(mass1d_min,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(delta1d_max,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(delta1d_min,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(z1d_max,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(z1d_min,1,MPI_DOUBLE_PRECISION,leader, tag1, MPI_COMM_WORLD, status, ierr)
		end if
	end if	!rank

end subroutine init_grizzly

!!This will read halofiles..
!!
subroutine part1(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
	integer :: ierr, i, tag1, tag3, status(MPI_STATUS_SIZE), haloid

	tag1=1
	tag3=3

	call det_leader(gpid, npgp, leader, gpext)

	if(rank ==leader) then
		comm_redshift=z_checkpoint(redshift_loop_index)
		write(*,*) ''
		write(*,*) 'Processing redshift, # = ',comm_redshift, redshift_loop_index
		call show_time('This slice starts at time: ')
		haloid=0
		if(genxhimap==1 .or. genalphamap==1) call read_halofile(rank, comm_redshift, num_halo, halomassnew, haloagenew, haloid)

		DO i=leader+1,gpext
		     	call MPI_SEND(num_halo, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
		     	call MPI_SEND(comm_redshift, 1, MPI_DOUBLE_PRECISION, i, tag1, MPI_COMM_WORLD, ierr) 
		END DO
	elseif(rank>leader .and. rank<=gpext) then
		call MPI_RECV(num_halo,1,MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(comm_redshift,1,MPI_DOUBLE_PRECISION, leader, tag1, MPI_COMM_WORLD, status, ierr)
	end if	!rank


end subroutine  part1

!!This reads the halo list, estimate ages of the sources. 
!!allocate  dat_overlap array which keep all the infor about the sources
!!
subroutine part2(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
	integer :: ierr, i, tag1, tag3, status(MPI_STATUS_SIZE), haloid

	call det_leader(gpid, npgp, leader, gpext)

	if(num_halo>=1) then
	if(rank==leader) then
		if(genxhimap==1) allocate(dat_overlap(6,num_halo))	!x, M, loss
		if(genxhimap==1) allocate(parameter1d_array(5,num_halo))	!!M1d, d1d, z1d, age1d, Vfinal 
		haloid=1
		matrix_nhi=1.0	!!this is to apply feedback correctly in the first step	!!this is required to prevent exit for selected redshofits
		if(genxhimap==1 .or. genalphamap==1) call read_halofile_track_age(rank, comm_redshift, num_halo, halomassnew, haloagenew, haloid)
	elseif(rank>leader .and. rank<=gpext) then
		if(genxhimap==1) allocate(dat_overlap(6,num_halo))
		if(genxhimap==1) allocate(parameter1d_array(5,num_halo))
	end if	!rank
	end if

end subroutine part2

!!This reads density field,
!!This determines the ionized bubble sizes and put them around the sources..
!!This step do not distribute the unused photons among the sources.
!!
subroutine part3(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
      	real(8) ::  mysum, mysum1, mysum2
	integer :: i, ierr, offset, onset, tag1, tag2, tag3, tag4, status(MPI_STATUS_SIZE), iid
	real(4), dimension(n_cell, n_cell, n_cell):: arrtmp

	tag1=1
	tag2=2
	tag3=3
	tag4=4
	call det_leader(gpid, npgp, leader, gpext)

!!First to read density field and share with all preocessor in the group

	if(num_halo>=1) then
	if(rank==leader) then
		call read_density(rank,comm_redshift)
		matrix_nhi=1.0	!!this is to apply feedback correctly in the first step
		matrix_tk=(tcmb*(1.0+comm_redshift)*(1.0+comm_redshift)/(z_dec))
		matrix_alpha= 0.0
		DO i=leader+1,gpext
			call MPI_SEND(over_density,n_cell*n_cell*n_cell,MPI_REAL,i,tag3, MPI_COMM_WORLD,ierr)
		END DO
		WRITE(*,*) 'Density field is read and shared with others.'
	elseif(rank>leader .and. rank<=gpext) then
		call MPI_RECV(over_density, n_cell*n_cell*n_cell, MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
	end if	!rank
	end if


	if(num_halo>=1) then
	if(rank==leader) then
		if(genxhimap==1) then
			DO i=leader+1,gpext
				iid=i-1-leader
				call para_range(num_halo, npgp-1, iid, offset, onset)	!!changes to exclude rank 0
			     	call MPI_SEND(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
			     	call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
				if(offset>0 .and. onset>=offset) then
			     		call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
			     		call MPI_SEND(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, ierr) 
				end if	!offset
			END DO
		end if

		if(genxhimap==1) then
			arrtmp=0.0

		  	DO i=leader+1,gpext
			     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
			     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
				if(offset>0 .and. onset>=offset) then
			     		call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
					arrtmp=arrtmp+matrix_nhi
				end if	!offset
			END DO

			call show_time('First step done: HII bubbles are created: redistribution remains')
			matrix_nhi=arrtmp

		  	DO i=leader+1,gpext
			     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
			     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
				if(offset>0 .and. onset>=offset) then
					call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag1, MPI_COMM_WORLD, ierr)
			     		call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
			     		call MPI_RECV(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, i, tag4, MPI_COMM_WORLD, status, ierr)
				end if	!offset
			END DO

			call mass_f(comm_redshift, global_expt_mf)

			if(global_expt_mf<0.05) matrix_nhi=0.0

			write(*,*) 'z, global_expt_mf', comm_redshift, global_expt_mf

			DO i=leader+1,gpext
			     	call MPI_SEND(global_expt_mf, 1, MPI_DOUBLE_PRECISION, i, tag1, MPI_COMM_WORLD, ierr) 
			END DO
		end if	!genxhimap
	elseif(rank>leader .and. rank<=gpext) then
		if(genxhimap==1) then
			call MPI_RECV(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, status, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
				call MPI_RECV(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, leader, tag4, MPI_COMM_WORLD, status, ierr)

				call cal_ion(offset,onset,rank,mysum1,comm_redshift)
				call visual_overlap1(offset,onset,rank,mysum2, matrix_nhi)
			end if

			call MPI_SEND(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, ierr) 
			call MPI_SEND(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
			end if


			call visual_overlap2(offset,onset,rank, matrix_nhi)

			call MPI_SEND(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, ierr) 
			call MPI_SEND(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag1, MPI_COMM_WORLD, status, ierr)
				call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
				call MPI_SEND(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, leader, tag4, MPI_COMM_WORLD, ierr)
			end if

			call MPI_RECV(global_expt_mf, 1, MPI_DOUBLE_PRECISION, leader, tag1, MPI_COMM_WORLD, status, ierr)
		end if	!genxhimap
	end if	!rank
	end if

end subroutine part3

!! This step to redistribute the unused photons among the sources.
!!
subroutine part4(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
      	real(8) ::  mf, redshift
	integer :: i, ierr, offset, onset, tag1, tag2, tag3, tag4, status(MPI_STATUS_SIZE)

	call det_leader(gpid, npgp, leader, gpext)

	if(num_halo>1) then
		tag1=1
		tag2=2
		tag3=3
		tag4=4

		if(visual_multicore) then
			call part4_0(gpid, npgp)
		else
			if(rank==leader) then
			if(genxhimap==1) then
				call visual_overlap3(1,num_halo,rank, matrix_nhi)
				parameter1d_array(5,1:num_halo) = dat_overlap(4,1:num_halo)
				call show_time('Second step done: redistribution of unused photons is done.')
				call mass_f(comm_redshift, mf)
			end if
			end if
		end if
	else
			if(rank==leader)	call show_time('Second step done: redistribution of unused photons is done.')
	end if

end subroutine part4




!!for visual 3
subroutine part4_0(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext


      	real(8) ::  mysum, mysum1, mysum2

	integer :: i, ierr, offset, onset, tag1, tag2, tag3, tag4, status(MPI_STATUS_SIZE), con_fstat, fstat, cn, ii,j,k, nn, iid
	real(4), dimension(n_cell, n_cell, n_cell):: arrtmp, matrix_nhi_back

	call det_leader(gpid, npgp, leader, gpext)

if(num_halo>=1 .and. global_expt_mf>0.05) then

	tag1=1
	tag2=2
	tag3=3
	tag4=4

	cn=0

DO ii=1,100

!!!overlap 1 2..
	if(rank==leader) then
	cn=cn+1

if(genxhimap==1) then


		matrix_nhi_back=matrix_nhi	!!xhii
		DO i=leader+1,gpext
			iid=i-1-leader
			call para_range(num_halo, npgp-1, iid, offset, onset)	!!changes to exclude rank 0

		     	call MPI_SEND(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
		     	call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
			if(offset>0 .and. onset>=offset) then
		     		call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
				call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr)
			end if	!offset
		END DO
end if



if(genxhimap==1) then
		arrtmp=0.0
		con_fstat=0

          	DO i=leader+1,gpext
		     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
			if(offset>0 .and. onset>=offset) then
		     		call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
		     		call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
		     		call MPI_RECV(fstat, 1, MPI_INTEGER, i, tag3, MPI_COMM_WORLD, status, ierr)
				arrtmp=arrtmp+matrix_nhi
				con_fstat=con_fstat+fstat
			end if	!offset

		END DO
 
		matrix_nhi=arrtmp

		DO i=1, n_cell
		DO j=1, n_cell
		DO k=1, n_cell
		if(matrix_nhi_back(i,j,k)>=1.0) matrix_nhi(i,j,k)=1.0

		END DO
		END DO
		END DO


          	DO i=leader+1,gpext
		call MPI_SEND(con_fstat, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
		END DO

	if(con_fstat==0) then
		parameter1d_array(5,1:num_halo) = dat_overlap(4,1:num_halo)
write(*,*) 'existing visual3 lop'
	 	exit
	end if

end if	!genxhimap


write(*,*) 'part4_0'

	elseif(rank>leader .and. rank<=gpext) then



if(genxhimap==1) then

		call MPI_RECV(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, status, ierr)
		if(offset>0 .and. onset>=offset) then
			call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
		     	call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
			call visual_overlap3_mod(offset,onset,rank, matrix_nhi, fstat)
		end if

		call MPI_SEND(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, ierr) 
		call MPI_SEND(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, ierr)
		if(offset>0 .and. onset>=offset) then
			call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
			call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
			call MPI_SEND(fstat, 1, MPI_INTEGER, leader, tag3, MPI_COMM_WORLD, ierr)
		end if

		call MPI_RECV(con_fstat, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)

		if(con_fstat==0) exit


end if	!genxhimap


	end if	!rank

END DO
end if

end subroutine part4_0


!!This account the Boundary effect of the HII bubble.
!!Grid pints will be partially ionized at the boundaries of the HII bubble.
!!
subroutine part5(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
      	real(8) ::  mysum, mysum1, mysum2, mf
	integer :: i, ierr, offset, onset, tag1, tag2, tag3, tag4, status(MPI_STATUS_SIZE), iid
	real(4), dimension(n_cell, n_cell, n_cell):: arrtmp

	tag1=1
	tag2=2
	tag3=3
	tag4=4

	call det_leader(gpid, npgp, leader, gpext)

	if(num_halo>=1 .and. global_expt_mf>0.05 .and. genxhimap==1) then
		if(rank==leader) then
			DO i=leader+1,gpext
				iid=i-1-leader
				call para_range(num_halo, npgp-1, iid, offset, onset)	!!changes to exclude rank 0
			     	call MPI_SEND(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, ierr) 
			     	call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
				if(offset>0 .and. onset>=offset) then
			     		call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
				end if	!offset
			END DO

			arrtmp=0.0
		  	DO i=leader+1,gpext
			     	call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
			     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
				if(offset>0 .and. onset>=offset) then
			     		call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
			     		call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, status, ierr)
					arrtmp=arrtmp+matrix_nhi
				end if	!offset
			END DO
			call bound_checkxhi(arrtmp)
			matrix_nhi=arrtmp
			call mass_f(comm_redshift,mf)
			
			call show_time('Third step: Accounted boundary effect.')
			WRITE(*,*) 'Max min xhii', maxval(matrix_nhi), minval(matrix_nhi)
!stop
		elseif(rank>leader .and. rank<=gpext) then
			call MPI_RECV(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, status, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_RECV(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
				call visual_overlap4(offset,onset,rank, matrix_nhi)
			end if

			call MPI_SEND(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, ierr) 
			call MPI_SEND(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
				call MPI_SEND(dat_overlap(1,offset),6*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, ierr)
			end if
		end if	!rank
	end if	!genxhimap


end subroutine part5

!!This will account the 1D profiles to the sources.
!!Also generate the TK map and lyalpha map
!!
subroutine part6(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
      	real(8) ::  mf, redshift
	integer :: i, ierr, offset, onset, tag1, tag2, tag3, tag4, status(MPI_STATUS_SIZE), iid, idps
	real(4), dimension(n_cell, n_cell, n_cell):: arrtmp, binaryxhi

	tag1=1
	tag2=2
	tag3=3
	tag4=4
	call det_leader(gpid, npgp, leader, gpext)

if(num_halo>=1 .and. global_expt_mf>0.05 .and. genxhimap==1) then
	if(rank==leader) then
		call mass_f(comm_redshift,mf)
		call xhii_xhi_binary(matrix_nhi, binaryxhi)
		call vol_f(comm_redshift,mf)
		call mass_f(comm_redshift,mf)


		if(useproilexhi==1) then
			call show_time('Fourth step: 1Dprofiles will be included')
			matrix_tk=matrix_nhi	!!tempraly saving xhi of the binary fields
			DO i=leader+1,gpext
				iid=i-1-leader
				call para_range(num_halo, npgp-1, iid, offset, onset)
				call MPI_SEND(offset, 1, MPI_INTEGER, i,  tag1, MPI_COMM_WORLD, ierr) 
				call MPI_SEND(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr)
				if(offset>0 .and. onset>=offset) then
					call MPI_SEND(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, i, tag3, MPI_COMM_WORLD, ierr) 
				end if
			END DO 
			call show_time('')

			if(genalphamap==1) call c_alpha_sim(rank, comm_redshift, redshift_loop_index, matrix_alpha)


			arrtmp=0.0
			DO i=leader+1,gpext
				call MPI_RECV(offset, 1, MPI_INTEGER, i, tag1, MPI_COMM_WORLD, status, ierr)
			     	call MPI_RECV(onset, 1, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, status, ierr)
				if(offset>0 .and. onset>=offset) then
			     		call MPI_RECV(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, i, tag4, MPI_COMM_WORLD, status, ierr)
					arrtmp=arrtmp + matrix_nhi*(1.0-arrtmp)	!!This step need to be corrected. need to incorprate already partially ionized fraction
				end if
			END DO
			binaryxhi=(1.0-arrtmp)*binaryxhi	!!xhi.. this do not include hii region at the boundary  
			matrix_nhi=binaryxhi*matrix_tk	!!this xhi array include hii region at the boundary 
		end if	! useproilexhi

		!!These steps are to produce TK, tb maps, PS etc

		if(gentkmap==1) then
			if(genxhimap==0) then
			 call read_xhifrac(rank, comm_redshift)
			binaryxhi = matrix_nhi
			end if	
				!write(*,*) 'May not produce correct effective temp map at the boundaries..'
			call gen_tkmap(rank, comm_redshift, binaryxhi, matrix_tk)
		end if	!gentkmap!

		redshift=z_checkpoint(redshift_loop_index+1)
		if(gentbmap==1) then
			if(genxhimap==0) call read_xhifrac(rank, redshift)
			if(genxhimap==0 .and. model_sc .ne. 0) call read_tkmatrix(rank, redshift)
			if(genxhimap==0 .and. model_sc == 3) call read_tkmatrix(rank, redshift)
			call delta_T(rank, redshift, over_density,  matrix_nhi, matrix_tk, matrix_alpha, arrtmp, model_sc)
		end if
		call write_data(redshift, matrix_nhi, matrix_tk, matrix_alpha, arrtmp)
		if(genps==1) then
			call psself(redshift, arrtmp,1)
		end if
		if(genpsxx==1)	then
			arrtmp=1.d0-matrix_nhi
			call psself(redshift, arrtmp,2)
		end if
		if(genpsxd==1)	then
			arrtmp=1.d0-matrix_nhi
			call pscross(redshift, arrtmp,over_density,1)
		end if
		call show_time('1Dprofiles are included')

	elseif(rank>leader .and. rank<=gpext) then
		if(useproilexhi==1) then
			call MPI_RECV(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
			call MPI_RECV(onset, 1, MPI_INTEGER,  leader, tag2, MPI_COMM_WORLD, status, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_RECV(parameter1d_array(1,offset),5*(onset-offset+1), MPI_REAL, leader, tag3, MPI_COMM_WORLD, status, ierr)
				call visual_nh_tk(offset,onset,rank, matrix_nhi)
			end if

			call MPI_SEND(offset, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, ierr) 
			call MPI_SEND(onset, 1, MPI_INTEGER, leader, tag2, MPI_COMM_WORLD, ierr)
			if(offset>0 .and. onset>=offset) then
				call MPI_SEND(matrix_nhi, n_cell*n_cell*n_cell, MPI_REAL, leader, tag4, MPI_COMM_WORLD, ierr)
			end if
		end if	!useproilexhi
	end if	!rank 
end if

end subroutine part6

!!This deallocate the allocated array
!!
subroutine part7(gpid, npgp)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer :: leader, gpext
      	real(8) ::  mf, redshift
	integer :: ierr, status(MPI_STATUS_SIZE)

	call det_leader(gpid, npgp, leader, gpext)
	if(num_halo>=1) then
		if(rank==leader) then
			if(genxhimap==1) deallocate(dat_overlap)
			if(genxhimap==1) deallocate(parameter1d_array)
		elseif(rank>leader .and. rank<=gpext) then
			if(genxhimap==1) deallocate(dat_overlap)
			if(genxhimap==1) deallocate(parameter1d_array)
		end if	!rank 
	end if
end subroutine part7

!!This checks for the end of the simulation
!!
subroutine part8(gpid, npgp, exit_id)
	use param
	implicit none
	include 'mpif.h'
	integer, intent(in) :: gpid, npgp
	integer, intent(inout) :: exit_id
	integer :: leader, gpext
      	real(8) ::  mf
	integer :: ierr, status(MPI_STATUS_SIZE), tag1, i


	tag1=1
	call det_leader(gpid, npgp, leader, gpext)

	if(num_halo>=1) then
		if(rank==leader) then
			if(parameter_space==0) then
				call vol_f(comm_redshift,mf)
				if(mf .gt. 0.990) exit_id=1
			end if

			DO i=leader+1,gpext
				call MPI_SEND(exit_id, 1, MPI_INTEGER, i,  tag1, MPI_COMM_WORLD, ierr) 
			END DO
		elseif(rank>leader .and. rank<=gpext) then
			call MPI_RECV(exit_id, 1, MPI_INTEGER, leader, tag1, MPI_COMM_WORLD, status, ierr)
		end if
	end if

	if(rank>=leader .and. rank<=gpext) then
		if(rank==leader)	call show_time('finish')
		if(exit_id==1 .and. useanalytic==0 .or. gentkmap==1) deallocate(ifront_data)
	end if
end subroutine part8


