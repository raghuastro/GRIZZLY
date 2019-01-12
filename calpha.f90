
!!This program calculate the x_alpha maps for a given source distribution
!!Assumption 1/r^2 fall of the Lyalpha photns number density: More details is in the README file
!!
subroutine c_alpha_sim(myid, z, zid, arral)
	use param
	implicit none 
	integer, intent(in) :: myid, zid
	real(8), intent(in) :: z
	real(4), dimension(n_cell, n_cell, n_cell), intent(out) :: arral
      	integer ::  i, j, k, check, cnt, z_loop_index_lya
      	real(8)::  af, sim_unit, zin
	real(4), dimension(n_cell, n_cell, n_cell):: dummy_arr



	if(stellar_sed) then
		call read_SED(nsed,lambdaarr, lumarr)
!	else
!		call generate_SED(nsed, lambdaarr, lumarr)
	end if


	check = zid

	lya_redshift=zid	!!This is the redshift where the coupling will be calculated..
	arral=0.0
	dummy_arr=0.0
	sim_unit=box/dble(n_cell)/hlittle/(1.0+lya_redshift) !mpc physical


!!start the loop to consider the cntributions from previous redshifts..
!DO z_loop_index_lya=check,check
DO z_loop_index_lya=1,check

		zin=z_checkpoint(z_loop_index_lya)
		!write(*,*) '***'
		!write(*,*) 'z loop start, index, z', z_loop_index_lya, zin
		gap_timestep=age_checkpoint(1)
		time_since= age_checkpoint(check)-age_checkpoint(z_loop_index_lya)+gap_timestep	!time gap between zstar and zstep
		imax_comm = int((time_since)*Megayear*c_light/Megaparsec/sim_unit)
		imin_comm = int((time_since-gap_timestep)*Megayear*c_light/Megaparsec/sim_unit)
		!write(*,*) 'timesince, timegap', time_since, gap_timestep
		!write(*,*) 'imin, imax, i, check', imin_comm, imax_comm, i, check



		call visual_calpha(zin, arral)


		!write(*,*) 'before:sum alpha this loop', sum(dummy_arr)
		dummy_arr=dummy_arr+arral
		!write(*,*) 'after:sum alpha this loop', sum(dummy_arr)


		!call cpu_time(time)
		!write(*,*) 'time ', time


	if(z_loop_index_lya==check) then
	
		arral=dummy_arr
		write(*,*) 'max, min, sum of alpha',maxval(arral), minval(arral), sum(arral)

	end if

end do

end subroutine c_alpha_sim




!!This include the SED for calculating the number of Lyamn alpha photons
!!SED is divided into two files. "lambda_mod.dat" contains all the lambda bins of the SED. Must contain the wavelengths below Lyman-alpha.
!! "lum_mod.dat" contains the Luminosity.
!! 'lambda' in angstrom, 'lum' in ev/s/lambda.. lambda should be in increasing order
!! One needs to replace this with prefered SED.. keeping the requirements ..
!!
SUBROUTINE read_SED(n_sed_e,lambda, lum)
	use param
	IMPLICIT NONE	
	integer, intent(in) :: n_sed_e
	real(8),dimension(n_sed_e), intent(out)::lambda,lum
  	character(180) :: ifile
	integer::i, fstat
	real(8)::xz, res, e , lambda_e
	integer, parameter :: n_sed_s=141		!!dim of the arry to read the SED

	if(n_sed_e .ne. n_sed_s) then
		write(*,*) 'Dimension of the input array is not matching, replace n_Sed to 141, input dimension', n_sed_e
		stop
	end if

	ifile=pspath//"lambda_mod.dat"
	!write(*,*) 'lambda array is read from ',ifile 
	fstat=0
	open(unit=21,file=ifile, status='old', iostat=fstat,  form='formatted')

	if(fstat .ne. 0) then
		write(*,*) 'File missing!!ERROR!!', ifile
		stop
	end if
	DO i=1,n_sed_s
		READ(21,*) xz
		lambda(i)=xz
	END DO
	close(21)

	ifile=pspath//"lum_mod.dat"
	!write(*,*) 'luminsity file array is read from ',ifile 
	open(unit=21,file=ifile, form='formatted', status='old', iostat=fstat)
	if(fstat /= 0) then
		write(*,*) 'File missing!!ERROR!!', ifile
		stop
	end if
	DO i=1,n_sed_s
		READ(21,*) xz
		lum(i)=xz
	END DO
	close(21)

!!normalised by the number of ionizing photons used in C2Ray
	res=0.0
	DO i=1, n_sed_s-1
		lambda_e=lambda(i)
		e=hplanck_ev*c_light/angstrom_cm/lambda_e	!!in ev
		if(e>E_HI) then
			res=res+lum(i)/e*(lambda(i+1)-lambda(i))*f_esc_hi
		end if
	END DO

	res=n_gamma_hi/res
	!write(*,*) 'CONVERSION FRACTOR TO MATCHE C2RAY IONIZING PHOTON RATE', res
	!write(*,*) 'SED Normalization is done to match the ionizing photons  per sec per solar mass halo (for g_gamma=1)', n_gamma_hi
	lum=lum*res


END SUBROUTINE read_SED




!! This subroutine calculates the Lyman-alpha coupling coefficient.
!! This takes the halos subset (myoffset,myonset) from rank 'myid'.
! 
subroutine visual_calpha(z, arral)
	use param
	use func
	IMPLICIT NONE
	real(8), intent(in):: z
	real(4), dimension(n_cell, n_cell, n_cell), intent(out) :: arral
	real(8):: r1, r2, al_min
	integer :: N, nh
	real(4), dimension(n_cell,n_cell,n_cell):: profilearr, arrhalo
	integer, dimension(n_cell, n_cell, n_cell) :: arrage
	!double precision, dimension(n_cell,n_cell,n_cell)::arrdbl
	double complex, dimension(n_cell,n_cell,n_cell) :: sourcearr_uv, profilearr_uv

	N=n_cell

	r1=dble(imin_comm)
	r2=dble(imax_comm)
if(r1>dble(n_cell/2)) then
	!write(*,*) 'Skipping this step as photons are out of box, z', z
	arral=0.0
else

	call read_halofile(0, z, nh, arrhalo, arrage, 1)
!write(*,*) 'max min M', maxval(arrhalo), minval(arrhalo)
	call generate_profile(n_cell, n_cell, n_cell, r1, r2, z, profilearr)
!write(*,*) 'max min profile', maxval(profilearr), minval(profilearr)

	!arrdbl=dble(profilearr)
	call FFTW_FORWARD_ARRAY(profilearr, profilearr_uv, N, N, N)
	!arrdbl=dble(arrhalo)
	call FFTW_FORWARD_ARRAY(arrhalo, sourcearr_uv, N,N, N)

	sourcearr_uv=sourcearr_uv*profilearr_uv

	call FFTW_BACKWARD_ARRAY(arral, sourcearr_uv, N, N, N)

	!arral=arr

	al_min=minval(arral)
	if(al_min<0.d0) arral=arral+al_min
!write(*,*) 'max min al', maxval(arral), minval(arral)
!stop

end if

end subroutine visual_calpha


!!This will create the 1/r/r profile for this particulater time step.. create profile between radius r1 and r2
!! r1 and r2 in grid unit..
!
SUBROUTINE generate_profile(n1, n2, n3, r1,r2,z, arr)
	IMPLICIT NONE
	integer, intent(in) ::n1, n2, n3
	real(8), intent(IN)::r1, r2, z
	real(4), dimension(n1, n2, n3), intent(out)::arr
	real(8)::ux, uy, uz, u, res, alpha
	integer::i,j,k

	DO i=1, n1
	DO j=1, n2
	DO k=1, n3
	ux=(dble(i-n1/2))+0.5d0
	uy=(dble(j-n2/2))+0.5d0
	uz=(dble(k-n3/2))+0.5d0
	u=sqrt(ux*ux+uy*uy+uz*uz)	!!grid unit
		if(U>r1 .and. U<r2) then
			call zeffectlyman_r(z,u,alpha)
			arr(i,j,k)=alpha
		else
			arr(i,j,k)=0.d0
		end if

	END DO
	END DO
	END DO

END SUBROUTINE  generate_profile

!! This subruitne will call the SED per solar mass and calculate the number of Lyman alpha photon 
!! This uses 1/r/r relation to calculate x_alpha
!!
SUBROUTINE zeffectlyman_r(z,rgrid,alpha)	!!this rgrid is in grid unit
	use param
	use func
	IMPLICIT NONE
	real(8),intent(in)::z, rgrid
	real(8),intent(out)::alpha
	real(8),dimension(nsed)::lambda,lum
	real(8)::r,zr,lr,l1,l0,e1,e0,en,nu_al,del_l,rp
	integer::i,j,k
	real(8)::xz,res,del_nu, nlya_emission


	lambda = lambdaarr
	lum = lumarr

	del_nu=c_light/angstrom_cm/lambda_alpha/lambda_alpha !!frequency difference for del lambda = 1 A
	nu_al=E_lyalpha/hplanck_ev !Lyman alpha freq in Hz

	alpha=0.d0


	r=rgrid*box/hlittle/dble(n_cell) !comoving distace at grid points in Mpc
	zr=rtoz(r,z)
	lr=lambda_alpha*(1.d0+zr)/(1.d0+z)	!!wavelength in the continuum SEd which redshifted to Lyman alpha at distance r
	if(lr<lambda_hi) then
		alpha = 0.d0	!!above EHI, photons will be absorbed for ionization
		!write(*,*) 'HI level reached at grid', i
	else 
			DO j=1,nsed
				if(lambda(j)>=lr) then
					l1=lambda(j)
					l0=lambda(j-1)
					e1=lum(j)
					e0=lum(j-1)
					exit
				end if
			END DO

			en=e0+ (lr-l0)*(e1-e0)/(l1-l0) !linear interpolation !energy per A per sec in erg
			del_l=(1.d0+zr)/(1.d0+z)! wave length differ for which  at r the wavelength differ will be 1 A!c_light/nu_al/nu_al*(1+z)/(1+zr)*1e8 !del lambda in A for unit frequency at r
			rp=r/(1.d0+zr)*megaparsec
			alpha=en*del_l/4.d0/pi/rp/rp/del_nu/E_lyalpha
			alpha=1.66d11/(1.d0+zr)*alpha*f_esc_lya*f_lya*(1.0+0.6*(1.0-f_esc_hi))  ! ly coeffi for unit mass

	end if

END SUBROUTINE zeffectlyman_r




