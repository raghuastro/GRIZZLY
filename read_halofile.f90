!!This subroutine reads the list of halo at redshift z

!!to track age, put track_age as true
!!The halo file should be given in grid positions, sum all haloes >Mcut at every grid points and place.
!!Sould not have many similar sources at the same grid points, otherwise age calculation will be wrong.
!!
subroutine read_halofile_track_age(myid, z, nhalo, halo_mass_new, halo_age_new)
	use param
	implicit none
	integer, intent(in)::myid
	real(8), intent(in)::z
	integer, intent(out) :: nhalo
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: halo_mass_new
	integer, dimension(n_cell, n_cell, n_cell), intent(inout):: halo_age_new
	integer :: id
	real(4), dimension(n_cell, n_cell, n_cell):: halo_mass_old
	integer, dimension(n_cell, n_cell, n_cell):: halo_age_old, halo_trace


!!choose the halo file : consider feedback, large/small source

	if(track_age==1) then
		halo_mass_old=halo_mass_new
		halo_age_old=halo_age_new
		write(*,*) 'max min halo mass (previous step)', maxval(halo_mass_new), minval(halo_mass_new)
		write(*,*) 'max min halo age (previous step)', maxval(halo_age_new), minval(halo_age_new)
	else
		halo_mass_old=0.0
		halo_age_old=0
	end if

	id=1
	call read_halofile(myid, z, nhalo, halo_mass_new, halo_age_new, id)

	if(track_age==1) then
		halo_trace=halo_age_new
		halo_age_new=halo_age_old+halo_age_new

		call check_halo_move(halo_mass_new, halo_age_new, halo_mass_old, halo_age_old, halo_trace)
	end if


!!Now we will calculate the effective mass of the halo over time, this is to account the mass evlution of the halo which is in general exponential.
!!However, we will not taken into account the dennisty contrast evolution around the source over the time and will work only with the current density distribution at that redshift.

	call prepare_halolist(nhalo, halo_mass_new, halo_age_new, halo_mass_old, halo_age_old)

	write(*,*) 'Halofile read succesfully.'
	write(*,*) 'max halo age (this step)', maxval(halo_age_new)

end subroutine read_halofile_track_age


subroutine prepare_halolist(nhalo, halo_mass_new, halo_age_new, halo_mass_old, halo_age_old)
	use param
	implicit none
	integer, intent(in) :: nhalo
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: halo_mass_new
	real(4), dimension(n_cell, n_cell, n_cell), intent(in) :: halo_mass_old
	integer, dimension(n_cell, n_cell, n_cell), intent(in):: halo_age_new, halo_age_old

	real(8) :: Meff
	integer:: i,j,k,cnt, xx


	dat_overlap=0.0	!!This array will carry the mass, X(:), etc.. 
	cnt=0


	DO k=1,n_cell
	DO j=1,n_Cell
	DO i=1,n_cell
		xx=halo_age_new(i,j,k)
		if(xx .ne. 0 ) then
			cnt=cnt+1
			Meff=(halo_mass_new(i,j,k)+halo_mass_old(i,j,k)*(real(halo_age_new(i,j,k))-1.))/real(halo_age_new(i,j,k))
			halo_mass_new(i,j,k)=Meff

			dat_overlap(1,cnt)=real(i)+0.5	
			dat_overlap(2,cnt)=real(j)+0.5	
			dat_overlap(3,cnt)=real(k)+0.5
			dat_overlap(4,cnt)=Meff 

if(track_age==0) then
			dat_overlap(5,cnt)=fix_age
else
			dat_overlap(5,cnt)=real(age_checkpoint(halo_age_new(i,j,k)))
end if

		end if
	END DO
	END DO
	END DO

end subroutine prepare_halolist


!put id =0 if only need to count nhalo, else put 1
!
subroutine read_halofile(myid, z, nhalo, halo_mass_new, halo_age_new, id)
	use param
	implicit none
	integer, intent(in)::myid, id
	real(8), intent(in)::z
	integer, intent(out) :: nhalo
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: halo_mass_new
	integer, dimension(n_cell, n_cell, n_cell), intent(inout):: halo_age_new


  	character(7) :: z_s
  	character(250) :: ifile
	real(8) :: halo_max, halo_min, M, Mass, mass_gal, M1
	integer::nhalo_max, i, j, k, ll, fstat, io

  	real(8), parameter :: Vol_m3 = (box/hlittle*(Megaparsec/100.0))**3.0
  	real(8), parameter :: M_grid = (omega_m*rho_ckgs*Vol_m3)/dble(nc**3.0)
  	real(8), parameter :: M_grid_s = M_grid/(M_sun/1.0d3) !mass unit in M_solar, particle mass=8*M_grid_s

!!choose the halo file : consider feedback, large/small source

	
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=halolist_path//trim(z_s)//trim(halofilename)
	open(unit=13,file=ifile, status='old')
	read(13, *) nhalo_max

	nhalo=0

	if(id .ne. 0) then
		halo_max=0.0
		halo_min=1.0d+20
		halo_mass_new=0.0		!!This is the new array of halo mass
		halo_age_new=0
	else
		write(*,*) 'halo file:', ifile
		write(*,*) 'Max number of halo', nhalo_max
	end if

!!Read the haloes :

	DO ll=1, nhalo_max ! halo loop
		read(13,*,IOSTAT=io)  i,j,k,M,M1
		if(io .ne. 0) exit
		!mass=M*M_grid_s*param_zeta   ! mass of the halo in solar mass
		mass=M*M_grid_s   ! mass of the halo in solar mass
		if(mass >=param_Mmin) then
			mass_gal=Mass*param_zeta
			nhalo=nhalo+1

		if(id .ne. 0) then
			if(mass_gal .lt. halo_min) halo_min=mass_gal
			if(mass_gal .gt. halo_max) halo_max=mass_gal

			halo_mass_new(i,j,k)=halo_mass_new(i,j,k)+real(mass_gal)
			halo_age_new(i,j,k)=1
		end if
		end if
				
	END DO !l
	!110 print*,'Halo considered',nhalo

	close(13)
	num_halo=nhalo

	if(id .ne. 0) then
		write(*,*) 'max min halo mass (this step)', halo_max, halo_min
	end if

end subroutine read_halofile



subroutine check_halo_move(halo_mass_new, halo_age_new, halo_mass_old, halo_age_old, halo_trace)
	use param
	implicit none
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: halo_mass_new, halo_mass_old
	integer, dimension(n_cell, n_cell, n_cell), intent(inout):: halo_age_new, halo_age_old
	integer, dimension(n_cell, n_cell, n_cell), intent(inout):: halo_trace

	integer:: i, j, k,cnt,nhalo_max, ll, xx, yy, cnt1, fstat
	integer::r1,r2,r3,n
	integer, dimension(3)::r

!!Check if halo is moving to different grid or not.. estimate the age of the halo..
	cnt=0
	cnt1=0
	DO k=1, n_cell
	DO j=1, n_cell
	Do i=1, n_cell
		xx=halo_age_old(i,j,k)
		yy=halo_age_new(i,j,k)

		if((xx .ne. 0) .and. (xx==yy) ) then

			!!Posibility of halo moving to new grid
				fstat=0
				DO n=2,100

					r(1)=i + cube(n,1)
					r(2)=j + cube(n,2)
					r(3)=k + cube(n,3)
					call check_boundary(r, n_cell, r1, r2, r3)

					if((halo_age_new(r1,r2,r3)==1) .and. (halo_trace(r1,r2,r3)==1) ) then !! possibility of halo move to new cell
						halo_age_new(r1,r2,r3)=halo_age_old(i,j,k) + 1	!!This halo will not be considered for next halo move, as the age have been changed to more than 1
						halo_mass_old(r1, r2, r3)=halo_mass_old(i,j,k)

						halo_mass_old(i,j,k)=0.0
						halo_age_new(i,j,k)=0
						halo_age_old(i,j,k)=0

						fstat=1
						exit
					end if
				END DO

			!!Possibility of halo merge to existing halo..

			if(fstat==0) then
				DO n=2,100
					r(1)=i + cube(n,1)
					r(2)=j + cube(n,2)
					r(3)=k + cube(n,3)
					call check_boundary(r, n_cell, r1, r2, r3)

					if( (halo_trace(r1,r2,r3)==1) ) then !! possibility of halo merge
					if(halo_age_new(r1,r2,r3) .lt. halo_age_new(i,j,k)) then
						halo_age_new(r1,r2,r3)=halo_age_new(i,j,k) + 1
						halo_age_old(r1,r2,r3)=halo_age_new(i,j,k)
					end if
						halo_age_new(i,j,k)=0
						halo_age_old(i,j,k)=0
						halo_mass_new(i,j,k)=0.0
						halo_mass_old(i,j,k)=0.0
						fstat=1
						exit
					end if
				END DO

			end if

			if(fstat==0) then
				halo_age_new(i,j,k)=0
				halo_age_old(i,j,k)=0
				halo_mass_new(i,j,k)=0.0
				halo_mass_old(i,j,k)=0.0
				fstat=1
			end if


			cnt1=cnt1+1
			if(fstat==0) cnt=cnt+1
		end if
	END DO
	END DO
	END DO


	!write(*,*) 'Number of halo moved, rejected', cnt1, cnt

end subroutine check_halo_move


