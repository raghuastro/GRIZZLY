SUBROUTINE show_time(ifile)
	IMPLICIT NONE
	character(*), intent(in) :: ifile
	real(4) :: time

	call cpu_time(time)
	write(*,*) ' Time:', time, ifile

END SUBROUTINE show_time


subroutine det_filename
use param
implicit none
	character(len=25)::  mmin_s
	character(len=7)::zeta_s, fx_s, al_s, sed_s, model_s,  age_s

	write(mmin_s,'(e25.3)') param_Mmin
	mmin_s=adjustl(mmin_s)
	write(zeta_s,'(f7.3)') param_zeta
	zeta_s=adjustl(zeta_s)
	write(age_s,'(I3)') track_age
	age_s=adjustl(age_s)

	write(sed_s,'(I3)') source_sed
	sed_s=adjustl(sed_s)
	write(model_s,'(I3)') model_sc
	model_s=adjustl(model_s)


	write(fx_s,'(f7.3)') param_fx
	fx_s=adjustl(fx_s)
	write(al_s,'(f7.3)') param_al
	al_s=adjustl(al_s)

	filename_comm="zeta"//trim(zeta_s)//"_Mmin"//trim(mmin_s)//"_sed"//trim(sed_s)//"_fx"//trim(fx_s) &
	& //"_al"//trim(al_s)//"_model"//trim(model_s)//"_ageT"//trim(age_s)

end subroutine det_filename


!!To describe the smulation..
!!
SUBROUTINE sim_specification
	use param
	IMPLICIT NONE

	if(parameter_space == 1) then
	if(mod(numtasks-1,gr_members_max) .ne. 0) then
	write(*,*) 'Choose number of processors smartly.'
	write(*,*) 'Number of group members', gr_members_max
	stop
	end if
	end if

	write(*,*) 'Cosmological parameters from WMAP5'
	write(*,*) 'Simulation box size (cMpc/h)=', box
	write(*,*) 'Simulation grid size=', n_cell
	write(*,*) 'Halo file will be read from..', halolist_path
	write(*,*) 'Output will be written in ..', output_path
	write(*,*) 'g_gamma for Large and small mass sources', g_gamma_c2ray_L, g_gamma_c2ray_S
	if(stellar_sed) then
	write(*,*) 'SED is choosen from PEGASE code.'
	else
	write(*,*) 'SED is choosen as blackbody'
	end if
	write(*,*) 'f_esc_hi=', f_esc_hi
	write(*,*) 'f_esc_lya=', f_esc_lya

	if(genxhimap==1 .and. writexhimap==1) write(*,*) 'Will generate xhi map and write data'
	if(genxhimap==1 .and. writexhimap==0) write(*,*) 'Will generate xhi map and not write data'
	if(genxhimap==0 .and. writexhimap==0) write(*,*) 'Will not generate xhi map and not write data'
	if(genxhimap==0 .and. writexhimap==1) then
	write(*,*) 'Will generate xhi map and write data'
	write(*,*) '!!ERROR!!choose right conditions'
	stop
	end if

	if(gentkmap==1 .and. writetkmap==1) write(*,*) 'Will generate tk map and write data'
	if(gentkmap==1 .and. writetkmap==0) write(*,*) 'Will generate tk map and not write data'
	if(gentkmap==0 .and. writetkmap==0) write(*,*) 'Will not generate tk map and not write data'
	if(gentkmap==0 .and. writetkmap==1) then
	write(*,*) 'Will generate xhi map and write data'
	write(*,*) '!!ERROR!!choose right conditions'
	stop
	end if

	if(genalphamap==1 .and. writealphamap==1) write(*,*) 'Will generate alpha map and write data'
	if(genalphamap==1 .and. writealphamap==0) write(*,*) 'Will generate alpha map and not write data'
	if(genalphamap==0 .and. writealphamap==0) write(*,*) 'Will not generate alpha map and not write data'
	if(genalphamap==0 .and. writealphamap==1) then
	write(*,*) 'Will generate xhi map and write data'
	write(*,*) '!!ERROR!!choose right conditions'
	stop
	end if

	if(gentbmap==1 .and. writetbmap==1) write(*,*) 'Will generate tb map and write data'
	if(gentbmap==1 .and. writetbmap==0) write(*,*) 'Will generate tb map and not write data'
	if(gentbmap==0 .and. writetbmap==0) write(*,*) 'Will not generate tb map and not write data'
	if(gentbmap==0 .and. writetbmap==1) then
	write(*,*) 'Will generate xhi map and write data'
	write(*,*) '!!ERROR!!choose right conditions'
	stop
	end if


	if(useanalytic==1 .and. useproilexhi == 1) then
	write(*,*) 'useanalytic==1 .and. useproilexhi == 1'
	write(*,*) '!!ERROR!!choose right conditions'
	stop
	end if

if(genps==1) then
	if(genalphamap==1 .and. model_sc .ne. 3) call error_msg('Alpha maps not required for model 1 or 2')
	if(model_sc==3 .and. genalphamap==0) call error_msg('Alpha maps required for model 3')
	if(model_sc==3 .and. gentkmap==0) call error_msg('Tk maps required for model 3')
	if(model_sc==2 .and. gentkmap==0) call error_msg('Tk maps required for model 2')
	if(gentkmap==1 .and. model_sc == 1) call error_msg('Tk maps not required for model 1')
end if



END SUBROUTINE sim_specification

SUBROUTINE error_msg(ifile)
IMPLICIT NONE
	character(*), intent(in) :: ifile
	write(*,*) ifile
	stop

END SUBROUTINE error_msg



SUBROUTINE mean_nh_nhe(z, nh, nhe)
	use param
	IMPLICIT NONE
	real(8)::z, nh, nhe

	nh=rho_ckgs*1d-3*omega_b/m_proton*hydrogen_frac*(1.+z)**3.0
	nhe=rho_ckgs*1d-3*omega_b/m_proton*(1.d0-hydrogen_frac)/4.d0*(1.+z)**3.0

	write(*,*) 'nh, nhe, z', nh, nhe, z

END SUBROUTINE mean_nh_nhe

!!This will optimize the distribution of jobs into many processes.
!
SUBROUTINE para_range(njob, nprocs, irank, myoffset, myonset)
	use param
	implicit none
	integer:: njob, nprocs, irank, myoffset, myonset,chunksize, i, nj, m

	if(njob<nprocs) then
		if(irank<njob) then
			myoffset=irank+1
			myonset=irank+1
		else
			myoffset=0
			myonset=0
		end if
	else
		myoffset=1
		nj=njob
		DO i=1, nprocs
			m=int(nj/(nprocs-i+1))
			myonset=myoffset+m-1
			if(irank==(i-1)) then
				exit
			else
			 	myoffset=myonset+1
				nj=nj-m
			end if
		END DO
	end if

	!write(*,*) 'JOB DISTRIBUTION, RANK, offset, onset', irank, myoffset, myonset

END SUBROUTINE para_range

!!This calculate the mass fraction of the ionization field
!
SUBROUTINE mass_f(z,mf)
	use param
	IMPLICIT NONE
	real(8),INTENT(IN)::z
	real(8),INTENT(OUT)::mf

	mf=1.d0 - sum(dble(MATRIX_NHI(:,:,:))*dble(over_density(:,:,:)))/(dble(n_cell*n_cell*n_cell))
	WRITE(*,*) 'Z=',z,'mass_frac=',mf

END SUBROUTINE mass_f

!!This calculate the volume fraction of the ionization field
!
SUBROUTINE vol_f(z,vf)
	use param
	IMPLICIT NONE
	real(8),INTENT(IN)::z
	real(8),INTENT(OUT)::vf

	vf=1.d0 - sum(dble(MATRIX_NHI(:,:,:)))/(dble(n_cell*n_cell*n_cell))
	WRITE(*,*) 'Z=',z,'vl_frac=',vf

END SUBROUTINE vol_f



subroutine cal_partial_ionfrac(z, vf)
use param
implicit none
real(8)::z, vf
integer::i,j,k, cnt1, cnt

	cnt=0
	DO i=1,n_cell
	DO j=1, n_cell
	DO k=1, n_cell
	if(matrix_nhi(i,j,k)>0.1) then
		cnt=cnt+1
	end if
	END DO
	END DO
	END DO
	 vf=dble(cnt)/dble(n_cell**3)

end subroutine cal_partial_ionfrac

subroutine hist_ionfrac(z, arr, xfracs)
use param
implicit none
	character(180)::ifile
	character(12)::z_s
real(8)::z, xmin, xmax,xhii
integer::i,j,k, ll
real(8), dimension(11, 2), intent(out)::arr
real(8), dimension(10), intent(out)::xfracs

arr=0.0
		write(z_s,'(f7.3)') z
		z_s=adjustl(z_s)
ifile=output_path//trim(z_s)//'histrogram_xhii.dat'
open(unit=190, file=ifile, status='replace', form='formatted')
DO i=1, size(arr(:,1))
arr(i,1)=dble(i-1)*0.1
END DO

DO i=1, n_cell
DO j=1, n_cell
DO k=1, n_cell
xhii=1.0-matrix_nhi(i,j,k)
DO ll=1, size(arr(:,1))-1
xmin=arr(ll,1)
xmax=arr(ll+1,1)
if(xhii>=xmin .and. xhii<xmax) then
arr(ll,2)=arr(ll,2) + 1.0
end if


END DO

END DO
END DO
END DO


DO ll=1, size(arr(:,1))-1
write(10,'(17f20.8)') arr(ll,1), arr(ll,2)/sum(arr(:,2))
xfracs(ll)=arr(ll,2)/sum(arr(:,2))
end do
	close(190)



end subroutine hist_ionfrac

!!This calculate the heated fraction of the ionization field
!
SUBROUTINE heated_frac(z,hf)
	use param
	IMPLICIT NONE
	real(8)::z,hf,tcmb_z
	integer::num,i,j,k

	tcmb_z=tcmb*(1.0+z)
	num=0
		DO k=1,n_cell
		DO j=1,n_cell
		DO i=1,n_cell
			if(matrix_tk(i,j,k) .gt. real(Tcmb_z)) then
			num=num+1
			end if
		END DO
		END DO
		END DO

	hf=dble(num)/dble(n_cell)**3.0
	!write(*,*) 'heated fraction=',hf
END SUBROUTINE heated_frac


!This write the three dimensional 4bit array to file in binary form
!
SUBROUTINE write_binary_4bit(arr, ifile, n1, n2, n3)
	IMPLICIT NONE
 	character(*)::ifile
	integer, intent(in)::n1, n2, n3
	real(4), dimension(n1, n2, n3)::arr
	integer::fstat

	open(unit=17,file=trim(ifile),status='replace',iostat=fstat,form='unformatted')
	!open(unit=17,file=trim(ifile),status='replace',iostat=fstat,form='binary')
	write(17) arr
	close(17)

	write(*,*) 'New file created :', ifile

END SUBROUTINE write_binary_4bit

!This reads 4bit three dimensional array from binary file
!
SUBROUTINE read_binary_4bit(arr, ifile)
	IMPLICIT NONE
 	character(*)::ifile
	real(4), dimension(:,:,:)::arr
	integer::fstat

	open(unit=17,file=trim(ifile),status='old',iostat=fstat,form='unformatted')
	!open(unit=17,file=trim(ifile),status='old',iostat=fstat,form='binary')
	read(17) arr
	close(17)

	write(*,*) 'New file created :', ifile

END SUBROUTINE read_binary_4bit

SUBROUTINE write_slice(arr, ifile, n1, n2, n3, n)
	IMPLICIT NONE
 	character(*)::ifile
	integer, intent(in)::n1, n2, n3, n
	real(4), dimension(n1, n2, n3)::arr
	integer::fstat,i,j

	real(4)::delr

	open(unit=17,file=trim(ifile),status='replace',iostat=fstat,form='formatted')

	DO i=1, n1
	DO j=1, n2
	write(17, '(17f20.4)') real(i), real(j), arr(i,j,n)
	END DO
	END DO

	 close(17)

END SUBROUTINE write_slice



!!z1, z2 are redshfit and delt is time gap in Myr
SUBROUTINE timegap(z1,z2, delt)
	use adaptint
	use param
	IMPLICIT NONE
	real(8)::z1,z2, delt

	double precision::RELACC,ABSACC,acc, res
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	interface
	FUNCTION fn_timegap(z)
	IMPLICIT NONE
	real(8),intent(in)::z
	real(8)::fn_timegap
	END FUNCTION fn_timegap
	end interface

	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0

	call D01ARF(z1,z2,fn_timegap,RELACC,ABSACC,MAXRUL,IPARM,ACC,res,N,&
		     &ALPHA,IFAIL)
	delt=res/megayear

END SUBROUTINE timegap


FUNCTION fn_timegap(z)
IMPLICIT NONE
	real(8),intent(in)::z
	real(8)::fn_timegap
	interface
	function hubble_constant(z) 
	use param
	implicit none
	real(8),intent(in)::z
	real(8)::hubble_constant
	end function hubble_constant
	end interface

	fn_timegap=1.d0/(1.d0+z)/hubble_constant(z) 


END FUNCTION fn_timegap

!!! return hubble costant at redshift z in sec-1
function hubble_constant(z) 
	use param
	implicit none
	real(8),intent(in)::z
	real(8)::hubble_constant
	hubble_constant=Ho*SQRT(omega_k*(1.d0+z)**2+omega_m*(1.d0+z)**3+&
         &omega_r*(1.d0+z)**4+omega_l)
end function hubble_constant


subroutine meminfo
	use param
	implicit none
	real(8)::mem
	integer::n, n4, n8

	n4=17 !! This is the number of 4 byte cube in the param file

	mem=dble(n_cell**3*n4) & 	!!These are the 4byte cubes like matrix_nhi..vel4 etc
	+ dble(max_cell_read*3) &	!!This is the 'cube' array
	+ dble(max_halos*12) &		!! This is the array which carry all information about the haloes.. this is the maxium array.. one should reduce the maximum halo number ot reduce the size of the array
	+ dble(n_ifront*5)		!!This is the ifront array 

	mem=mem*4.0	!!byte
	mem=mem/1024.0	!!KB
	mem=mem/1024.0	!!MB
	mem=mem/1024.0 	!!GB

	write(*,*) 'Memory required for the simulation : (GB)', mem*1.1 	
!stop
end subroutine meminfo 

SUBROUTINE pearson_cross_correlation(arr1, arr2, n1, n2, corel)
	IMPLICIT NONE
	integer, intent(in)::n1, n2
	real(8), dimension(n1, n2), intent(in)::arr1, arr2
	real(8), intent(out)::corel

	integer::i,j
	real(8)::av1, av2, x1, x2, x3


	av1=sum(arr1)/dble(n1*n2)
	av2=sum(arr2)/dble(n1*n2)
	x1=0.d0
	x2=0.d0
	x3=0.d0
	DO i=1, n1
	DO j=1, n2
		x1=x1 + (arr1(i,j)-av1)*(arr2(i,j)-av2)
		x2=x2 + (arr1(i,j)-av1)**2.d0
		x3=x3 + (arr2(i,j)-av2)**2.d0
	END DO
	END DO

	corel=x1/sqrt(x2*x3)

END SUBROUTINE pearson_cross_correlation



subroutine delta_T(myid, z, arrdd,  arrxhi, arrtk, arral, arrtb, model)
	use param
	implicit none
	integer, intent(in) :: myid, model
	real(8), intent(in) :: z
	real(4), dimension(n_cell, n_cell, n_cell), intent(in) :: arrdd, arrxhi, arrtk, arral
	real(4), dimension(n_cell, n_cell, n_cell), intent(out) :: arrtb


	real(8):: Tcmb_z


	!write(*,*) 'delta_T is called : rank ', myid

	Tcmb_z=tcmb*(1.0+z)
	arrtb=0.0

!!calculate spin temp 

	if(model==1) then
		arrtb=1e4
	elseif(model==2) then
		arrtb=arrtk
	elseif(model==3) then
		arrtb=arrtk
 		arrtb(:,:,:)=(1.0+arral(:,:,:))/((1.0/real(Tcmb_z))+(arral(:,:,:)/arrtk(:,:,:)))
	end if



	arrtb(:,:,:)=26.0*arrxhi(:,:,:)*real((omega_b*hlittle*hlittle/0.0230))* &
		real(sqrt((1.0+z)/10.0*(0.30/omega_m)))*(arrdd(:,:,:))*(1.0-real(Tcmb_z)/arrtb(:,:,:))

	!write(*,*) 'RANK,max min tb',myid,maxval(arrtb),minval(arrtb)


end subroutine delta_T


subroutine write_data(z, arrxhi, arrtk, arral, arrtb)
	use param
	implicit none
	real(8), intent(in) :: z
	real(4), dimension(n_cell, n_cell, n_cell), intent(in) :: arrxhi, arrtk, arral, arrtb
  	character(180) ::  ofile
  	character(7) :: z_s
		write(z_s,'(f7.3)') z
		z_s=adjustl(z_s)


			if(writexhimap==1) then
				ofile=output_path//trim(z_s)//trim(filename_comm)//trim(xhifilename)
				call write_binary_4bit(arrxhi, ofile, n_cell, n_cell, n_cell)
			end if	!genxhimap

			if(writetkmap==1) then
				ofile=output_path//trim(z_s)//trim(filename_comm)//trim(tkfilename)
				call write_binary_4bit(arrtk, ofile, n_cell, n_cell, n_cell)
			end if	!gentkmap

			if(writealphamap==1) then
				ofile=output_path//trim(z_s)//trim(filename_comm)//trim(alphafilename)
				call write_binary_4bit(arral, ofile, n_cell, n_cell, n_cell)
			end if	!gentkmap
			if(writetbmap==1) then
				ofile=output_path//trim(z_s)//trim(filename_comm)//trim(tbfilename)
				call write_binary_4bit(arrtb, ofile, n_cell, n_cell, n_cell)
			end if


end subroutine write_data



subroutine psself(z, arr, id)
	use param
	use, intrinsic :: iso_c_binding
	IMPLICIT NONE
	include 'fftw3.f03'
	real(8), INTENT(IN)::z
	real(4), DIMENSION(n_cell, n_cell, n_cell), INTENT(IN)::arr
	integer, intent(in) :: id	!1 Tb, 2 x

	INTEGER, dimension(Nbins_ps)::Nptsarr
	real(8), dimension(Nbins_ps):: pskarr, kbinarr

	INTEGER::myid, N1, N2, N3, Nbins
	REAL(8):: boxsize

	real(8), DIMENSION(n_cell, n_cell, n_cell)::MAP

	INTEGER ::i, j, k, ii, jj, kk, ll, l,  fstat
	REAL(8)::k_min,k_max,k_x,k_y,k_z,k_vec, volume_box, pixelsize, mean, avtb, k3, rms, skewness, kurtosis, dtb, global_mean_tb, mf
	character(180) :: ofile



	REAL(8)::pk, delk1, delk2, grid3d, grid2d, cons, delk
	character(len=7) z_s
	real(8),dimension(n_cell,n_cell)::slice

	integer*8 p, q
	double complex, dimension(n_cell, n_cell, n_cell) :: out_data, in_data



	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)
	if(id==1) then
		ofile=output_path//trim(z_s)//trim(filename_comm)//trim(pskfilename)
	elseif(id==2) then
		ofile=output_path//trim(z_s)//trim(filename_comm)//trim(pskxxfilename)
	else
		write(*,*) 'chose proper ps scenario:terminating program'
		stop
	end if

	MAP=dble(arr)

	Nbins=Nbins_ps
	N1=n_cell
	N2=n_cell
	N3=n_cell


	boxsize=box/hlittle

	call init_ps_log(Nbins_ps, kbinarr)


	grid3d=dble(N1*N2*N3)
	volume_box=(boxsize)**3.0
	pixelsize=volume_box/grid3d

	k_min=kbinarr(1)
	k_max=kbinarr(Nbins)


	call dfftw_plan_dft_3d(p, N1, N2, N3, in_data, out_data, FFTW_FORWARD, FFTW_ESTIMATE)



	global_mean_tb=sum(MAP)/grid3d
	grid2d=dble(N1*N2)
	DO k=1, N3
		mean=sum(MAP(:,:, k))/grid2d !! 2D mean
		slice(:,:)=MAP(:,:,k)-mean
		MAP(:,:,k)=slice(:,:)
		in_data(:,:,k)=dcmplx(slice(:,:),0.0)
	END DO

	!write(*,*) 'Mean is subtracted from each channel'


	call dfftw_execute_dft(p, in_data, out_data)

	pskarr=0.0
	Nptsarr=0

	cons=2.0*pi/boxsize
	DO i=1,N1
	DO j=1,N2
	DO k=1,N3
		ii=i
		jj=j
		kk=k
		if(i>N1/2) ii=i-N1
	 	if(j>N2/2) jj=j-N2
	 	if(k>N3/2) kk=k-N3

		k_x=dble(ii-1)*cons
		k_y=dble(jj-1)*cons
		k_z=dble(kk-1)*cons
		k_vec=sqrt(k_x*k_x+k_y*k_y+k_z*k_z)

		pk=(dreal(out_data(i,j,k))*dreal(out_data(i,j,k)) + dimag(out_data(i,j,k))*dimag(out_data(i,j,k)))*pixelsize*pixelsize/volume_box ! //power spectrum unit Mpc^3


		DO l=1,Nbins
			if(l==1) then
				delk1=(kbinarr(2)-kbinarr(1))*0.5d0
				delk2=delk1
			elseif(l==Nbins) then
				delk1=(kbinarr(Nbins)-kbinarr(Nbins))*0.5d0
				delk2=delk1
			else
				delk1=(kbinarr(l)-kbinarr(l-1))*0.5d0
				delk2=(kbinarr(l+1)-kbinarr(l))*0.5d0
			end if

			if(k_vec>=kbinarr(l)-delk1/2.0 .AND. k_vec<kbinarr(l)+delk2/2.0) then
					Nptsarr(l)=Nptsarr(l)+1
					pskarr(l)=pskarr(l)+pk
			end if

		END DO

	END DO
	END DO
	END DO

		call dfftw_destroy_plan(p)



		DO l=1,Nbins_ps
			k_vec=kbinarr(l)
			k3=k_vec**3.d0/2.d0/pi/pi
			if(Nptsarr(l)>0) pskarr(l)=pskarr(l)/dble(Nptsarr(l))*k3
		END DO


!!other quantities of interest

	mean=sum(MAP)/dble(n_cell**3)
	rms=0.d0
	skewness=0.d0
	kurtosis=0.d0
	DO i=1, N1
	DO j=1, N2
	DO k=1, N3
	dtb=MAP(i,j,k)-mean
	rms=rms+dtb**2.d0
	skewness=skewness+dtb**3.d0
	kurtosis=kurtosis+dtb**4.d0
	END DO
	END DO
	END DO
	rms=rms/dble(N1*N2*n3)
	rms=sqrt(rms)
	skewness=skewness/dble(N1*N2*N3)/rms**3.d0
	kurtosis=kurtosis/dble(N1*N2*N3)/rms**4.d0


	call vol_f(z,mf)

	write(*,*) 'PSk file', ofile, global_mean_tb, rms, skewness, kurtosis
	open(unit=1115,file=trim(ofile),status='replace',iostat=fstat,form='formatted')
	write(1115,*) '#global mean of the quantity, xhii, rms, skewness, kurtosis'
	write(1115,*) global_mean_tb, mf, rms, skewness, kurtosis
	write(1115,*) '#k(1/Mpc), Del^2'
	DO l=1,Nbins_ps

		if(Nptsarr(l)>0) WRITE(1115,'(17f20.4)') kbinarr(l), pskarr(l), dble(Nptsarr(l))
		flush(1115)
	END DO
	close(1115)



!stop

end subroutine psself






subroutine pscross(z, arr1, arr2, id)
	use param
	use, intrinsic :: iso_c_binding
	IMPLICIT NONE
	include 'fftw3.f03'
	real(8), INTENT(IN)::z
	real(4), DIMENSION(n_cell, n_cell, n_cell), INTENT(IN)::arr1, arr2
	integer, intent(in) :: id

	INTEGER, dimension(Nbins_ps)::Nptsarr
	real(8), dimension(Nbins_ps):: pskarr, kbinarr

	INTEGER::myid, N1, N2, N3, Nbins
	REAL(8):: boxsize

	real(8), DIMENSION(n_cell, n_cell, n_cell)::MAP

	INTEGER ::i, j, k, ii, jj, kk, ll, l,  fstat
	REAL(8)::k_min,k_max,k_x,k_y,k_z,k_vec, volume_box, pixelsize, mean, avtb, k3, rms, skewness, kurtosis, dtb, global_mean_tb, mf
	character(180) :: ofile



	REAL(8)::pk, delk1, delk2, grid3d, grid2d, cons, delk, mean1, mean2
	character(len=7) z_s
	real(8),dimension(n_cell,n_cell)::slice


	integer*8 p, q
	double complex, dimension(n_cell, n_cell, n_cell) :: p_out_data, q_out_data,  p_in_data, q_in_data



	Nbins=Nbins_ps
	N1=n_cell
	N2=n_cell
	N3=n_cell


	boxsize=box/hlittle

	call init_ps_log(Nbins_ps, kbinarr)

	grid3d=dble(N1*N2*N3)
	volume_box=(boxsize)**3.0
	pixelsize=volume_box/grid3d

	k_min=kbinarr(1)
	k_max=kbinarr(Nbins)


	call dfftw_plan_dft_3d(p,N1,N2,N3, p_in_data, p_out_data, FFTW_FORWARD, FFTW_ESTIMATE)

	!for first cube
	MAP=dble(arr1)
	mean1=sum(map(:,:,:))/grid3d

	grid2d=dble(N2*N3)
	DO i=1, N1
		mean=sum(map(i,:,:))/grid2d !! 2D mean
		slice(:,:)=map(i,:,:)-mean
		p_in_data(i,:,:)=dcmplx(slice(:,:),0.0)
	END DO

	call dfftw_execute_dft(p, p_in_data, p_out_data)


	call dfftw_plan_dft_3d(q,N1,N2,N3, q_in_data, q_out_data, FFTW_FORWARD, FFTW_ESTIMATE)
	


!second cube
	MAP=dble(arr2)
	mean2=sum(map(:,:,:))/grid3d

	DO i=1, N1
		mean=sum(map(i,:,:))/grid2d !! 2D mean
		slice(:,:)=map(i,:,:)-mean
		q_in_data(i,:,:)=dcmplx(slice(:,:),0.0)
	END DO
	

	call dfftw_execute_dft(q, q_in_data, q_out_data)

	 

	pskarr=0.0
	Nptsarr=0

	cons=2.0*pi/boxsize
	DO i=1,N1
	DO j=1,N2
	DO k=1,N3
		ii=i
		jj=j
		kk=k
		if(i>N1/2) ii=i-N1
	 	if(j>N2/2) jj=j-N2
	 	if(k>N3/2) kk=k-N3

	!the factor 2 pi scales the power spectrum for exp(ik r) convention from  exp(i2 pi k r) convention, see notes for details*/

		k_x=dble(ii-1)*cons
		k_y=dble(jj-1)*cons
		k_z=dble(kk-1)*cons
		k_vec=sqrt(k_x*k_x+k_y*k_y+k_z*k_z)

		pk=(dreal(p_out_data(i,j,k))*dreal(q_out_data(i,j,k)) + dimag(p_out_data(i,j,k))*dimag(q_out_data(i,j,k)))*pixelsize*pixelsize/volume_box ! //power spectrum unit Mpc^3

		DO l=1,Nbins
			if(l==1) then
				delk1=(kbinarr(2)-kbinarr(1))*0.5d0
				delk2=delk1
			elseif(l==Nbins) then
				delk1=(kbinarr(Nbins)-kbinarr(Nbins))*0.5d0
				delk2=delk1
			else
				delk1=(kbinarr(l)-kbinarr(l-1))*0.5d0
				delk2=(kbinarr(l+1)-kbinarr(l))*0.5d0
			end if

			if(k_vec>=kbinarr(l)-delk1/2.0 .AND. k_vec<kbinarr(l)+delk2/2.0) then
					Nptsarr(l)=Nptsarr(l)+1
					pskarr(l)=pskarr(l)+pk
			end if

		END DO


	END DO
	END DO
	END DO


		call dfftw_destroy_plan(p)
		call dfftw_destroy_plan(q)

		DO l=1,Nbins_ps
			k_vec=kbinarr(l)
			k3=k_vec**3.d0/2.d0/pi/pi
			if(Nptsarr(l)>=1) pskarr(l)=pskarr(l)/dble(Nptsarr(l))*k3
		END DO


	write(z_s,'(f7.3)') z
	z_s=adjustl(z_s)
	if(id==1) then
		ofile=output_path//trim(z_s)//trim(filename_comm)//trim(pskxdfilename)
	else
		stop
	end if
	write(*,*) ofile
	
	open(unit=15,file=ofile,status='replace',iostat=fstat,form='formatted')

	DO l=1,Nbins_ps
		if(Nptsarr(l)>=1) WRITE(15,'(17f20.4)') kbinarr(l), pskarr(l), dble(Nptsarr(l))
		flush(1115)
	END DO

 	close(15)





end subroutine pscross




subroutine linspace_dp(n,xmin,xmax,x)
    	implicit none
	integer, intent(in) :: n
    	real(8),intent(in) :: xmin,xmax
    	real(8),dimension(n),intent(out) :: x
    	integer :: i

    	if (n == 1) then
	       	if(xmin /= xmax) then
		  	write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
		  	stop
	       	else
        		x = xmin
       		end if
    	else
       		do i=1,n
          		x(i) = (xmax-xmin) * dble(i) / dble(n-1) + xmin
       		end do
    	end if
end subroutine linspace_dp

subroutine logspace_dp(n,xmin,xmax,x)
    	implicit none
	integer, intent(in) :: n
    	real(8),intent(in) :: xmin,xmax
    	real(8), dimension(n),intent(out) :: x
	real(8) :: a, b
	double precision :: aa
	integer :: i
    	real(8), dimension(n) ::y
   	 if (size(x) == 1 .and. xmin /= xmax) then
       		write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       		stop
    		end if
    	call linspace_dp(n,log10(xmin),log10(xmax),y)

DO i=1, n
	a=y(i)

	aa=10.d0**dble(abs(a))
if(a>=0.0) then
    	x(i) = aa
else
    	x(i) = 1.0/aa
end if
END DO

end subroutine logspace_dp


!This will divide the k range into equal gaps
!
subroutine init_ps(Nbin, kbin)
	use param
	implicit none
	integer, intent(in) :: Nbin
	real(8), dimension(Nbin) :: kbin	
	integer::l
	real(8)::k_min,k_max,delk,boxsize

	boxsize=box/hlittle
	k_min=2.0*pi/boxsize*sqrt(3.0) !Mpc^-1
	k_max=k_min*(n_cell/2.0)*sqrt(3.0) ! sqrt(3.) comes because k=sqrt(kx^2+ky^2+kz^2)

	WRITE(*,*) 'K_min=',k_min
	WRITE(*,*) 'k_max=',k_max 

	call linspace_dp(Nbins_ps,k_min,k_max,kbin)

	!write(*,*) 'kbin max,min',maxval(kbin),minval(kbin)

end subroutine init_ps

!!This will divide k range into log binning
!
subroutine init_ps_log(Nbin, kbin)
	use param
	implicit none
	integer, intent(in) :: Nbin
	real(8), dimension(Nbin) :: kbin
	integer::l
	real(8)::k_min,k_max,delk,boxsize
	real(8), parameter::delkfix=1.5

	boxsize=box/hlittle
	k_min=2.d0*pi/boxsize*sqrt(3.0) !Mpc^-1
	k_max=k_min*(n_cell/2.0) ! sqrt(3.) comes because k=sqrt(kx^2+ky^2+kz^2)

	WRITE(*,*) 'K_min=',k_min
	WRITE(*,*) 'k_max=',k_max 

	call logspace_dp(Nbins_ps,k_min,k_max,kbin)

	!write(*,*) 'kbin max,min',maxval(kbin),minval(kbin)
	!write(*,*) 'kmin is increased'

end subroutine init_ps_log


SUBROUTINE write_slice4bit(arr, ifile, n1, n2, n3, n)
IMPLICIT NONE
 	character(*)::ifile
	integer, intent(in)::n1, n2, n3, n
	real(4), dimension(n1, n2, n3)::arr
	integer::fstat,i,j

	real(4)::delr

	write(*,*) ifile
	open(unit=17,file=trim(ifile),status='replace',iostat=fstat,form='formatted')

	DO i=1, n1
	DO j=1, n2
	write(17, '(17f20.4)') real(i), real(j), arr(i,j,n)
	END DO
	END DO

	 close(17)


END SUBROUTINE write_slice4bit


subroutine FFTW_FORWARD_ARRAY(MAP, MAPuv, N1, N2, N3)
	use param
	use, intrinsic :: iso_c_binding
	IMPLICIT NONE
	include 'fftw3.f03'
	INTEGER, INTENT(IN):: N1, N2, N3
	real(4), DIMENSION(N1, N2, N3), INTENT(IN)::MAP
	DOUBLE COMPLEX, DIMENSION(N1, N2, N3), INTENT(out):: MAPuv

	integer::p
	type(C_PTR) :: p_idata, p_odata
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: out_data
	complex(C_DOUBLE), dimension(:,:,:), pointer :: in_data


	p_idata = fftw_alloc_complex(int(N1 * N2 * N3,  C_SIZE_T))
	p_odata = fftw_alloc_complex(int(N1 * N2 * N3,  C_SIZE_T))

	call c_f_pointer(p_idata, in_data, [N1,N2,N3])
	call c_f_pointer(p_odata, out_data, [N1,N2,N3])

	call dfftw_plan_dft_3d(p, N1, N2, N3, in_data, out_data, FFTW_FORWARD, FFTW_ESTIMATE)

	in_data(:,:,:)=dcmplx(MAP,0.0)

	call dfftw_execute_dft(p, in_data, out_data)

	MAPuv=out_data


	call dfftw_destroy_plan(p)
	call fftw_free(p_idata)
	call fftw_free(p_odata)

end subroutine FFTW_FORWARD_ARRAY


SUBROUTINE FFTW_BACKWARD_ARRAY(MAP, MAPuv, N1, N2, N3)
	use, intrinsic :: iso_c_binding
	IMPLICIT NONE
	include 'fftw3.f03'
	integer, intent(in)::N1, N2, N3
	real(4), dimension(N1, N2, N3), intent(out):: MAP
	DOUBLE complex, dimension(N1, N2, N3), intent(in):: MAPuv

	integer::p
	type(C_PTR) :: p_idata, p_odata
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: out_data
	complex(C_DOUBLE), dimension(:,:,:), pointer :: in_data



	p_idata = fftw_alloc_complex(int(N1 * N2 * N3,  C_SIZE_T))
	p_odata = fftw_alloc_complex(int(N1 * N2 * N3,  C_SIZE_T))

	call c_f_pointer(p_idata, in_data, [N1,N2,N3])
	call c_f_pointer(p_odata, out_data, [N1,N2,N3])


	map=0.0
	call dfftw_plan_dft_3d(p,N1,N2,N3, out_data, in_data, FFTW_BACKWARD, FFTW_ESTIMATE)
	out_data=dcmplx(MAPuv)
	call dfftw_execute_dft(p, out_data, in_data)
	MAP=dreal(in_data)/dble(N1*N2*N3)

	call dfftw_destroy_plan(p)
	call fftw_free(p_idata)
	call fftw_free(p_odata)

END SUBROUTINE FFTW_BACKWARD_ARRAY






