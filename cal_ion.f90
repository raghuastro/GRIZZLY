

!!Calculate the ionized region, overdensity around a source::: given mass of the source is known

SUBROUTINE cal_ion(myoffset,myonset,myid,mysum1,z)
	use param
	implicit none
	integer ::myoffset,myonset,myid,i,l,fstat
	real(8),dimension(3)::halo_cm
	real(8) :: mass_gal,mysum1,av_delta,z_c,m_c,d_c,gal_life,ifront,trial_radius,z,sum_t, t_c, vol
  	character(180) :: ifile,ofile
  	character(7) :: z_s
	integer :: imin, imax

!!need to be updated to make it faster... z, fx, al, M, del, age, ifront
	!write(*,*) 'Calling cal_ion..',myoffset,myonset,myid,mysum1,z
call det_ifrontmin_imax(z, imin, imax)
!imin=1
!imax=n_ifront
!write(*,*) imin, imax, z
!stop

	mysum1=0.0

if((myoffset >0) .and. (myonset>=myoffset)) then
	DO i=myoffset,myonset
!if(myoffset==1) write(*,*) i, myoffset, myonset, myid, num_halo
		halo_cm(1:3)=dble(dat_overlap(1:3,i))
		mass_gal=dble(dat_overlap(4,i))
		gal_life=dble(dat_overlap(5,i))

		trial_radius=0.50*box/hlittle/(1.0+z)/dble(n_cell)

		DO l=1,500  
			call calculate_av_density(halo_cm,trial_radius,z,av_delta) 
			if(useanalytic==0) call ifront_mp(z,mass_gal,av_delta,ifront,d_c,m_c,z_c,t_c,gal_life, imin, imax)
			if(useanalytic==1) call stromgen(z,mass_gal,av_delta,ifront,d_c,m_c,z_c,t_c,gal_life)

			if(ABS((ifront-trial_radius)/ifront) .lt. 1e-3) then
				exit
			end if
			if(ifront .gt. trial_radius ) then
				trial_radius=trial_radius+(ifront-trial_radius)/2.0
			end if

			if(ifront .lt. trial_radius ) then
				trial_radius=ifront/1.2
			end if

		END DO !l

		parameter1d_array(1,i)=real(m_c) 
		parameter1d_array(2,i)=real(d_c)
		parameter1d_array(3,i)=real(z_c)
		parameter1d_array(4,i)=real(t_c)

		vol=trial_radius*(1.0+z)*dble(n_cell)*hlittle/box
		vol=(4.0/3.0)*pi*vol**3.d0

		dat_overlap(4,i)=real(vol) !!ionized volume around the source

		
		mysum1=mysum1+vol


	END DO

end if

	mysum1=mysum1/(dble(n_cell))**3.0
	! WRITE(*,*) 'called cal_ion.id,vol frac=',myid,mysum1

END SUBROUTINE cal_ion

SUBROUTINE det_ifrontmin_imax(z, imin, imax)
use param
IMPLICIT NONE
real(8), intent(in) ::z
integer, intent(out) :: imin, imax
real(8) :: zp, z1, z2
integer :: i, fstat

imin=1
imax=n_ifront
z1=100.0
z2=100.0
fstat=0
	DO i=1,n_ifront

		zp=dble(ifront_data(i,1))


		if((z2.ne.zp) .and. (zp<z)) then
z2=zp
fstat=fstat+1

if(fstat==2) then
		imax=i
		exit
end if
		end if

		if((zp .ne. z1) .and. (z<zp)) then
		z1=zp
		imin=i
		end if

	END DO

!write(*,*) z, dble(ifront_data(imin,1)), dble(ifront_data(imax,1))
!stop


end SUBROUTINE det_ifrontmin_imax

!!Calculate the average density around a halo with center 'halo_cm in sim unit', within a radius 'trail_radius in physical mpc', av_delta --> gm/cm3/<av rho>
!!
subroutine calculate_av_density(halo_cm,trial_radius,redshift,av_delta)
	use param
	implicit none
	real(8),dimension(3), intent(in)::halo_cm  !in global coordinate
	real(8)::trial_radius,redshift,av_delta  !in mpc
	integer :: count_cell1,q

	real(8)::vol1,sim_unit_length,max_cell,unit_vol
	integer,dimension(3)::r,x
	integer::r1,r2,r3

	vol1=0.0
	count_cell1=0
	sim_unit_length=box/dble(n_cell)/(1.0+redshift)/hlittle !mpc
	unit_vol=sim_unit_length**3.0

	x(:)=int(halo_cm(:)) 

	max_cell=(4.0*pi/3.0*trial_radius**3.0)/unit_vol
	if(max_cell < 1.0) then
		max_cell=1.0
	end if

	DO q=1,nint(max_cell)
		r(:)=x(:)+cube(q,:)
		call check_boundary(r, n_cell, r1, r2, r3)
		vol1=vol1+dble(over_density(r1,r2,r3)) 
		count_cell1=count_cell1+1
	END DO

	av_delta=vol1/dble(count_cell1)

end subroutine calculate_av_density

SUBROUTINE stromgen(z,mass,deltain,ifront,d_c,m_c,z_c,t_c,gal_life)
	use param
	IMPLICIT NONE
	real(8),intent(in)::z
	real(8),intent(in)::mass  !in solar mass
	real(8),intent(in)::deltain ! in gm/cm^3
	real(8),intent(in)::gal_life
	real(8),intent(out)::ifront
	real(8),intent(out)::d_c,m_c,z_c, t_c
	real(8)::n_gamma_dot, T, Alpha_B, co_clumpling, n_hidrogen, delta, t_rec, r_s, t_time, r_1, m_star, age, tot_phot, tot_nh
	integer::i


	delta=deltain!gas_density/av_rho_unnorm 

	n_gamma_dot=1.33d43*mass

!!Calculate the recombination time scale..

	T=1d4

	Alpha_B=2.59d-13*(T/1d4)**(-0.7d0)
	co_clumpling=1.0!2.9d0*((1.d0+z)/6.d0)**(-1.1d0)!1.d0  
	n_hidrogen=delta*((1.0d0+z)**3.d0)*1.87d-7
	t_rec=1.d0/(co_clumpling*Alpha_B*n_hidrogen)/megayear !in sec

!!calculate the stromgen radius using the ionizationn rate from the source..

	r_s=(3.d0*n_gamma_dot/(4.d0*pi*Alpha_B*co_clumpling*n_hidrogen*n_hidrogen))**(1.d0/3.d0)/kiloparsec	!!1d3 for kiloparsec

	t_time=gal_life
	r_1=r_s*(1.d0-exp(-(t_time/t_rec)))**(1.d0/3.d0)

	ifront=r_1/1d3

tot_phot=n_gamma_dot*gal_life*megayear
tot_nh=n_hidrogen*4.0/3.0*pi*(ifront*megaparsec)**3.0

if(tot_phot<tot_nh) then
write(*,*) 'Oh no', tot_phot, tot_nh, z,mass,deltain
stop
end if


m_star=mass
!delta=deltain
age=gal_life


	if(m_star < mass1d_min) call msg_outbound(m_star, mass1d_min, 'ERROR : mass is less than')
	if(m_star > mass1d_max) call msg_outbound(m_star, mass1d_max, 'ERROR : mass is larger than')
	if(delta < delta1d_min) call msg_outbound(delta, delta1d_min, 'ERROR : delta is less than')
	if(delta > delta1d_max) call msg_outbound(delta, delta1d_max, 'ERROR : delta is larger than')
	if(age < age1d_min) call msg_outbound(age, age1d_min, 'ERROR : age is less than')
	if(age > age1d_max) call msg_outbound(age, age1d_max, 'ERROR : age is larger than')


	



	d_c=dble(nint(delta))
	i=int(log10(m_star))
	m_c=10.**(dble(i))	!!need to correct
	z_c=dble(int(z))
	t_c=dble(nint(age/10.0)*10)



END SUBROUTINE stromgen

SUBROUTINE msg_outbound(M, Mb, ifile)
	IMPLICIT NONE
	real(8), intent(inout) :: M, Mb
 	character(*) :: ifile

	M=Mb 


END SUBROUTINE msg_outbound




!! This subroutine takes the mass of the qso and the constant density around it and 
!! create the ionised radius in Mpc within the QSO lifetime using photon conservation.
!! mass in solar mass, gas_density
subroutine ifront_mp(z,mass,deltain,ifront,d_c,m_c,z_c,t_c,gal_life, imin, imax)
	use param
	implicit none

	real(8),intent(in)::z, mass, deltain, gal_life
	real(8),intent(out)::ifront, d_c,m_c,z_c,t_c
	integer, intent(in) :: imin, imax


	real(8)::delta,den_c,m_star,d, ifront1, I0, I1, I2, age, zp
	real(8)::a1,a2,a3,a4,t1
	integer::fstat,i,count, con

	con=n_ifront
	delta=deltain
	m_star=mass
	age=gal_life
	zp=z

	if(m_star < mass1d_min) call msg_outbound(m_star, mass1d_min, 'ERROR : mass is less than')
	if(m_star > mass1d_max) call msg_outbound(m_star, mass1d_max, 'ERROR : mass is larger than')
	if(delta < delta1d_min) call msg_outbound(delta, delta1d_min, 'ERROR : delta is less than')
	if(delta > delta1d_max) call msg_outbound(delta, delta1d_max, 'ERROR : delta is larger than')
	if(age < age1d_min) call msg_outbound(age, age1d_min, 'ERROR : age is less than')
	if(age > age1d_max) call msg_outbound(age, age1d_max, 'ERROR : age is larger than')
	if(zp < z1d_min) call msg_outbound(zp, z1d_min, 'ERROR : z is less than')
	if(zp > z1d_max) call msg_outbound(zp, z1d_max, 'ERROR : z is larger than')


	DO i=imin,imax!1,n_ifront
		a1=dble(ifront_data(i,5))
		a2=dble(ifront_data(i,6))
		a3=dble(ifront_data(i,1))
		t1=dble(ifront_data(i,4))
		a4=dble(ifront_data(i,7))

		if(a1 >= m_star) then
		if(a2 >=delta) then	
		if(a3 >= zp) then
		if(t1>=age) then
			I0=a4/1000.0*((1.0+a3)/(1.0+z)) *(mass/a1)**(1.0/3.0) *(a2/deltain)**(1.0/3.0)*(gal_life/t1)**(1./3.) 
			d_c=a2
			m_c=a1
			z_c=a3
			t_c=t1
			exit
		end if
		end if
		end if
		end if






		if(i==n_ifront) then
			write(*,*) 'ERROR!!end of file'
			write(*,*) 'z,mass,deltain,ifront,d_c,m_c,z_c,gal_life', z,mass,deltain,ifront,d_c,m_c,z_c,gal_life
			write(*,*) 'from ifront file. m, d, z, t', a1, a2, a3, t1
			write(*,*) 'm_star, nint(d), nint(z)', m_star, nint(delta), nint(z)
			!stop
		end if
	END DO

	DO i=1,n_ifront
		a1=dble(ifront_data(i,5))
		a2=dble(ifront_data(i,6))
		a3=dble(ifront_data(i,1))
		t1=dble(ifront_data(i,4))
		a4=dble(ifront_data(i,7))


		if(a1 >= m_star) then
		if(a2 >= delta) then	
		if(a3 < zp) then
		if(t1>=age) then
			I1=a4/1000.0*((1.0+a3)/(1.0+z)) *(mass/a1)**(1.0/3.0) *(a2/deltain)**(1.0/3.0)*(gal_life/t1)**(1./3.) 
			exit
		end if
		end if
		end if
		end if





		if(i==n_ifront) then
			write(*,*) 'ERROR!!end of file'
			write(*,*) 'z,mass,deltain,ifront,d_c,m_c,z_c,gal_life', z,mass,deltain,ifront,d_c,m_c,z_c,gal_life
			write(*,*) 'from ifront file. m, d, z, t', a1, a2, a3, t1
			write(*,*) 'm_star, nint(d), nint(z)', m_star, nint(delta), nint(z)
			stop
		end if
	END DO


	ifront=(I1 + (I0-I1)*(z-real(int(z))) + I0-(I0-I1)*(real(int(z)+1)-z))/2.0

!if(m_c<1d9) then
!	ifront=ifront*0.9
!end if!

end subroutine ifront_mp

!!Read ifront file in array : mass, delta, z,time,ifront
!!This will be used for determining the 1D profile
!!
SUBROUTINE read_ifront(myid)
	use param
	implicit none
	integer::myid,con,fstat,i
	real(8), dimension(nparam1d)::arr
	character(len=180) ifile, ifrontfile, cha
	character(7) sed_s

	write(sed_s,'(I3)') source_sed
	sed_s=adjustl(sed_s)

	ifrontfile=output_1d//'./ifront_data_sed'//trim(sed_s)//'.dat'
	write(*,*) 'Ifront data will be read from',  ifrontfile

  	ifile=ifrontfile
  	open(unit=41,file=ifile,status='old',iostat=fstat,form='formatted')
	read(41,*) cha
	read(41,*) con
	n_ifront=con

!write(*,*) n_ifront, nparam1d, myid
!stop
	allocate(ifront_data(n_ifront, nparam1d))
  	open(unit=41,file=ifile,status='old',iostat=fstat,form='formatted')
	DO i=1,con
		read(41,'(17f20.8)') arr
		ifront_data(i,:)=real(arr)
	END DO
  	close(41)



	z1d_max=maxval(ifront_data(1:con,1))
	z1d_min=minval(ifront_data(1:con,1))
	fx1d_max=maxval(ifront_data(1:con,2))
	fx1d_min=minval(ifront_data(1:con,2))
	al1d_max=maxval(ifront_data(1:con,3))
	al1d_min=minval(ifront_data(1:con,3))
	age1d_max=maxval(ifront_data(1:con,4))
	age1d_min=minval(ifront_data(1:con,4))
	mass1d_max=maxval(ifront_data(1:con,5))
	mass1d_min=minval(ifront_data(1:con,5))
	delta1d_max=maxval(ifront_data(1:con,6))
	delta1d_min=minval(ifront_data(1:con,6))



	write(*,*) 'Ifront is read. Max min mass=',mass1d_max,mass1d_min
	write(*,*) 'Ifront is read. Max min delta=',delta1d_max,delta1d_min
	write(*,*) 'Ifront is read. Max min z=',z1d_max,z1d_min
	write(*,*) 'Ifront is read. Max min time=',age1d_max,age1d_min
	write(*,*) 'Ifront is read. Max min fx=',fx1d_max,fx1d_min
	write(*,*) 'Ifront is read. Max min alpha=',al1d_max,al1d_min
	write(*,*) 'Ifront is read. Max min ifront=',maxval(ifront_data(:,nparam1d)),minval(ifront_data(:,nparam1d))

END SUBROUTINE read_ifront


