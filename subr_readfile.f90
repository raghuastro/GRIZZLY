!!!!read the halofinds redshifts
!
subroutine read_redshift
	use param
	implicit none
	integer::fstat,i
	real(4)::arg

	open(11,file=checkpoints,status='old',iostat=fstat)
	read(11,*) i
	num_checkpoints=i

	DO i=1,num_checkpoints
	    read(11,*) arg
		z_checkpoint(i)=dble(arg)
	END DO

if(selected_reion ==0) then
	write(*,*) 'halofinds to recompose:'
	  do i=1,num_checkpoints
	    write(*,'(f5.1)') z_checkpoint(i)
	  enddo
end if


end subroutine read_redshift

subroutine read_select
	use param
	implicit none
	integer::fstat,i
	real(4)::arg

	open(11,file=selectpoints,status='old',iostat=fstat)
	read(11,*) num_selectpoints

	DO i=1,num_selectpoints
	    read(11,*) arg
		z_selectpoint(i)=dble(arg)
	END DO


	write(*,*) 'selected z to recompose:'
	  do i=1,num_selectpoints
	    write(*,'(f5.1)') z_selectpoint(i)
	  enddo


end subroutine read_select

subroutine cal_age_timestep
	use param
	implicit none
	real(8)::z1, z2, delt
	integer::i

	z2=z_checkpoint(1)

	if(num_checkpoints==1) then
		age_checkpoint(num_checkpoints)=10.0
	else

		DO i=1,num_checkpoints-1
		z1=z_checkpoint(i+1)
		call timegap(z1,z2,delt)
		age_checkpoint(i)=delt
		END DO
		age_checkpoint(num_checkpoints)=age_checkpoint(num_checkpoints-1) + 10.0
	end if

	!write(*,*) 'Age to recompose:'
	!do i=1,num_checkpoints
	!    	write(*,'(f5.1)') age_checkpoint(i)
	!enddo

end subroutine cal_age_timestep


!!This will generate an array of nearby grid points from 0,0,0 within a box of dimension n,n,n
!!
SUBROUTINE readcube_new(n)
	use param
	implicit none
	integer, intent(in) ::n
	integer :: nin
	real(4), dimension(:), allocatable::arr, arr1
	integer, dimension(:,:), allocatable ::cube1
	integer::i,j,k,ii,jj,kk,cnt, stat, cn, fstat
	real(4) ::d, tmp
	integer, dimension(3)::temp
	character(180) :: ifile
	character(12) :: grid_s

	nin=n/2
	allocate(arr(nin*nin*nin), arr1(nin*nin*nin), cube1(nin*nin*nin,3))

	ii=1
	jj=1
	kk=1
	cnt=0

	DO i=1, nin
	DO j=1, nin
	DO k=1, nin
		cnt=cnt+1
		d=real(i-ii)**2.0 + real(j-jj)**2.0 + real(k-kk)**2.0
		d=sqrt(d)
		arr(cnt) = d

		cube1(cnt,1) = i-ii
		cube1(cnt,2) = j-jj
		cube1(cnt,3) = k-kk
	END DO
	END DO
	END DO

	arr1=arr


	!write(*,*) 'Before'
	!write(*,*) arr(1), arr(2), maxval(arr)

	call sort(nin*nin*nin,arr)

	!write(*,*) 'After'
	!write(*,*) arr(1), arr(2), arr(cnt)

	cn=0
	DO i=1, cnt
		!if(mod(i,10000)==0) write(*,*) i
		d=arr(i)
		temp=cube1(i,:)
		tmp=arr1(i)
		stat=0
		DO j=i,cnt
		if(d==arr1(j)) then
			ii=cube1(j,1)
			jj=cube1(j,2)
			kk=cube1(j,3)
			stat=1
			arr1(j) =1e30
			arr1(i)=arr1(j)
			arr1(j)=tmp
			cube1(j,:)=temp


			cn=cn+1
			cube(cn,:)=(/-ii,-jj,-kk/)

			if(kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/-ii,-jj,kk/)
			end if
			if(jj .ne. -jj) then
				cn=cn+1
				cube(cn,:)=(/-ii,jj,-kk/)
			end if
			if(ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,-jj,-kk/)
			end if
			if(jj .ne.  -jj .and. kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/-ii,jj,kk/)
			end if

			if(ii .ne.  -ii .and. kk .ne. -kk) then
				cn=cn+1
				cube(cn,:)=(/ii,-jj,kk/)
			end if
			if(jj .ne.  -jj .and. ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,jj,-kk/)
			end if

			if(jj .ne.  -jj .and. kk .ne. -kk .and. ii .ne. -ii) then
				cn=cn+1
				cube(cn,:)=(/ii,jj,kk/)
			end if

			exit
		end if
		END DO
		if(stat==0) then
		write(*,*) 'not convergeing'
		stop
		end if
	END DO

	write(grid_s,'(I3)') ncube
	grid_s=adjustl(grid_s)


	ifile=cellpath//'cube'//trim(grid_s)//'.bi'
	write(*,*) 'cube array writting to file', ifile
	open(unit=311,file=ifile, form='unformatted', status='replace', iostat=fstat)
	write(311) cn
	write(311) cube
	close(311) 
	spdim=cn
	write(*,*) 'cube array dim', cn, n*n*n, n, real(cn)**(1./3.)

	write(*,*) cube(1,:), cube(2,:), cube(cn,:)
	deallocate(arr, arr1, cube1)

END SUBROUTINE readcube_new


!! This reads the array of nearby points from 0,0,0.. in case file misiing will generate the array
!!
SUBROUTINE readcube
	use param
	implicit none
  	character(180) ::ofile5,ifile1, ifile
	integer::q,q1,count,i,j,k,fs,n2p,count1,num, fstat
	character(12) :: grid_s

	write(grid_s,'(I3)') ncube
	grid_s=adjustl(grid_s)


	ifile=cellpath//'cube'//trim(grid_s)//'.bi'
	open(unit=31,file=ifile, form='unformatted', status='old', iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!..will be regenerated', ifile
		call	readcube_new(ncube)
	else
		read(31) q
		spdim=q
		if(q .ne. max_cell_read) then
		write(*,*) 'Dimension of cube array not matching.. terminating job.', q, max_cell_read
		stop
		end if
		read(31) cube
	end if	

	close(31)

	!write(*,*) 'done cube reading, maxval cube, spdim',maxval(cube), cube(1,:), spdim
!stop

END SUBROUTINE readcube


subroutine read_density(myid,z)
	use param
	implicit none
  	character(7) :: z_s
	integer::fstat,myid
	real(8)::z
  	character(180) ::ifile
	integer::n1, n2, n3
	real(8) :: av_rho_unnorm, norm_density

	REAL(8), parameter :: den_unit0=omega_m*rho_ckgs*(dble(n_cell)/dble(nc))**3.0	!Density 


	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=denisty_path//trim(z_s)//trim(densityfilename)
	open(unit=12,file=ifile,status='old',iostat=fstat,form='binary')
	read(12) n1, n2, n3
	read(12) over_density !dm density in sim unit
	close(12) 

	norm_density= den_unit0*(omega_b/omega_m) *(1.0+z)*(1.0+z)*(1.0+z)*1e-3 ! in gm/cm^3
	av_rho_unnorm=dble(sum(over_density(:,:,:)))/(dble(n_cell))**3.0



	write(*,*) 'Dimension of density cube', n1, n2, n3
	WRITE(*,*) '<rho>(gmcm-3)  sim, cosmology =',av_rho_unnorm*norm_density,(omega_b)*rho_ckgs*1d-3*(1.0+z)**3.0

	over_density=over_density/real(av_rho_unnorm)
	av_rho_unnorm=1.0

end subroutine read_density

subroutine read_xhifrac(myid,z)
	use param
	implicit none
  	character(7) :: z_s
	integer::fstat,myid
	real(8)::z
  	character(180) ::ifile
	integer :: n1, n2, n3

	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)


	ifile=output_path//trim(z_s)//trim(filename_comm)//trim(xhifilename)
write(*,*) ifile
	open(unit=31,file=ifile,status='old',iostat=fstat,form='unformatted')
	read(31) matrix_nhI
	close(31)

	write(*,*) 'max min xhi', maxval(matrix_nhi), minval(matrix_nhi)

end subroutine read_xhifrac

subroutine read_tkmatrix(myid,z)
	use param
	implicit none
  	character(7) :: z_s
	integer::fstat,myid, c2rayid
	real(8)::z, tkmin
  	character(180) ::ifile

	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=output_path//trim(z_s)//trim(filename_comm)//trim(tkfilename)
	open(unit=31,file=ifile,status='old',iostat=fstat,form='unformatted')
	read(31) matrix_tk
	close(31)

	write(*,*) 'max min TK', maxval(matrix_tk), minval(matrix_tk)
	tkmin=2.73*(1.+z)**2.0/(z_dec)
	write(*,*) 'Min Tk must be larger than', tkmin

end subroutine read_tkmatrix


subroutine read_tbmatrix(myid,z, arr)
	use param
	implicit none
	real(4), dimension(n_cell, n_cell, n_cell), intent(out) :: arr
  	character(7) :: z_s
	integer::fstat,myid, c2rayid
	real(8)::z, tkmin
  	character(180) ::ifile

	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=output_path//trim(z_s)//trim(filename_comm)//trim(tbfilename)
	write(*,*) ifile
	open(unit=31,file=ifile,status='old',iostat=fstat,form='unformatted')
	read(31) arr
	close(31)

	write(*,*) 'max min Tb', maxval(arr), minval(arr)


end subroutine read_tbmatrix



subroutine read_alphamatrix(myid,z)
	use param
	implicit none
  	character(7) :: z_s
	integer::fstat,myid, c2rayid
	real(8)::z
  	character(180) ::ifile

	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=output_path//trim(z_s)//trim(filename_comm)//trim(alphafilename)
	open(unit=31,file=ifile,status='old',iostat=fstat,form='unformatted')
	read(31) matrix_alpha
	close(31)

	write(*,*) 'max min x_alpha', maxval(matrix_alpha), minval(matrix_alpha)

end subroutine read_alphamatrix

subroutine check_boundary(x, n, i, j, k)
	implicit none
	integer, dimension(3), intent(in)::x
	integer, intent(in)::n
	integer, intent(out)::i,j,k

	i=x(1)
	j=x(2)
	k=x(3)

	if(i<1) i=i+n
	if(j<1) j=j+n
	if(k<1) k=k+n

	if(i>n) i=i-n
	if(j>n) j=j-n
	if(k>n) k=k-n

end subroutine check_boundary

subroutine read_velocity(myid,z, arr)
	use param
	implicit none
	real(4), dimension(3,n_cell,n_cell,n_cell), intent(out)::arr
  	character(7) :: z_s
	integer::fstat,myid, n1, n2, n3, n4
	real(8)::z, len_unit, tau_t, vel_unit
  	character(180) ::ifile
	real(4), dimension(n_cell, n_cell, n_cell) :: rho_gas


	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)

	ifile=denisty_path//z_s(1:len_trim(z_s))//'v_all.dat'
	open(unit=12,file=ifile,status='old',iostat=fstat,form='unformatted')
	read(12) n1, n2, n3
	write(*,*) n1, n2, n3
	read(12) arr
	close(12) 

	ifile=denisty_path//z_s(1:len_trim(z_s))//'n_all.dat'
	open(unit=12,file=ifile,status='old',iostat=fstat,form='binary')
	read(12) n1, n2, n3
	read(12) rho_gas !dm density in sim unit
	close(12) 


	arr(1,:,:,:)=arr(1,:,:,:)/rho_gas(:,:,:)*8.0
	arr(2,:,:,:)=arr(2,:,:,:)/rho_gas(:,:,:)*8.0
	arr(3,:,:,:)=arr(3,:,:,:)/rho_gas(:,:,:)*8.0

	len_unit=box*megaparsec/hlittle/(1.0+z)/dble(nc)
	tau_t=2.0/3.0/sqrt(omega_m*Ho*Ho)/(1.0+z)/(1.0+z)
	vel_unit=len_unit/tau_t

	arr=arr*vel_unit	!!cm/s rho/8

	write(*,*) 'max min velocity', maxval(arr), minval(arr)

end subroutine read_velocity




!! This read a 1d xhi and tk profile around a source and use that to make a correlation 
!! to generate the TK and xhi maps

!!This also read the ionization maps generated by overlap_mpi

SUBROUTINE input_rt_file(m, d, z, t, sed, fx, alpha, filename)
	use param
	IMPLICIT NONE
	real(8), intent(in)::m, d, z, fx, alpha, t
	integer, intent(in)::sed
 	character(*), intent(out)::filename
  	character(7) :: z_s, d_s, yr_s, sed_s, fx_s, al_s
  	character(25) :: m_s

	write(m_s,'(e25.3)') m
	m_s=adjustl(m_s)
	write(d_s,'(f7.3)') d
	d_s=adjustl(d_s)
	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)
	write(sed_s,'(I3)') sed
	sed_s=adjustl(sed_s)
	write(fx_s,'(f7.3)') fx
	fx_s=adjustl(fx_s)
	write(al_s,'(f7.3)') alpha
	al_s=adjustl(al_s)

	write(yr_s,'(I3)') int(t)
	yr_s=adjustl(yr_s)

	filename=output_1d//"m"//trim(m_s)//"_d"//trim(d_s)//"_sed"//trim(sed_s)//"_fx"//trim(fx_s) &
	& //"_al"//trim(al_s)//"_z"//trim(z_s)//'_t'//trim(yr_s)//'.bin'


END SUBROUTINE input_rt_file

subroutine read_reion(myid, z)
	use param
	implicit none
  	character(7) :: z_s
  	character(180) :: ifile, ofile
	integer::i, fstat, myid, ii, jj, k
	real(8)::T0, Tg, z, xhii


	logical,parameter::write_corel=.false.
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	real(4), dimension(:,:), allocatable::arr_bi
	real(8) :: M, d, t, red, Tk

	real(8), dimension(n_corel) :: xhiibin
	real(8), parameter :: xhimin = 1d-4
	real(8), parameter :: xhimax = 1.d0

	!write(*,*) 'read_reion is called : rank,z',myid,z


	M= mass1d_max
	d=delta1d_min
	t=age1d_max
	red=dble(int(z))!z1d_min

	call logspace_dp(n_corel,xhimin,xhimax,xhiibin)

	DO i=1, n_corel
		corel(i,1) = xhiibin(n_corel-i+1)
	END DO


!! read corel
	call input_rt_file(M, d, red, t, source_sed, param_fx_fix_1d, param_al_fix_1d, ifile)
	!call input_rt_file(1d10, 1.d0, dble(int(z)), 200.d0, source_sed, 0.05d0, 1.5d0, ifile)
	write(*,*) 'Corelation is read from file:', ifile
	!open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='binary')
	open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='unformatted')

	if (fstat /= 0) write(*,*) 'error opening file',ofile
		
	read(14) ii,jj
	allocate(arr_bi(ii,jj))
	read(14) arr_bi
	close(14)


	Tg=tcmb*(1.0+z)
	T0=Tg*(1.0+z)/(z_dec)
	corel(:,2)=real(T0)

	k=1
	DO i=1,ii  !grid
		xhii=1.d0-arr_bi(i, 2)
		if(xhii<=corel(k,1)) then
			Tk=arr_bi(i, 3)
			if(Tk<T0) Tk=T0
			corel(k,2) = Tk
			k=k+1
		end if
		if(k>n_corel) exit
		if(Tk<=T0) exit
	END DO
	


	if(write_corel) then
		ofile=output_path//'corel_file.dat'
		open(unit=315,file=ofile)
		write(*,*) 'Corrlation is written in file:',ofile

		DO i=1,n_corel
			write(315,*) corel(i,:)
		END DO
		close(315)
		stop
	end if

	deallocate(arr_bi)

	!write(*,*) 'Rethink about the binning used in corel.. may be log bin'

end subroutine read_reion


