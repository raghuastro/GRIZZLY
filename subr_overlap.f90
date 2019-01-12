

!! This one will generate the HII bubbles with overlap. photon not conserved.
!! This also creat partial ionization fraction of the HII bubble is very small.
!!
SUBROUTINE visual_overlap4(myoffset, myonset, myid, arr)
	use param
	implicit none
	integer, intent(in) :: myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr
	INTEGER:: q,q2,r1,r2,r3,nn,cells,con, n1, n2, fstat, den, i,j,k
	INTEGER,DIMENSION(3)::x, r
	real(8), dimension(3)::rr
	real(8):: extra,mysum,mf, rad, v, radhii, xhii, delr, ra, mfexp, ra1

	con=0
	extra=0.0
	mf=0.0
	nn=n_cell
	arr=0.0
!write(*,*) 'Overlap 4 is called', myoffset, myonset, myid, num_halo

!if(myid==1) write(*,*) 'Overlap 4 is called!!!THIS CAN BE IMPROVED BY COMING FROM THE BACK OF CUBE ARRAY KEPING INPUT MATRIXNHI' 

if((myoffset >0) .and. (myonset>=myoffset)) then

	DO q=myoffset,myonset
		x(1:3)=int(dat_overlap(1:3,q)) 
		radhii=(3.0/4.0/pi*dat_overlap(4,q))**(1.0/3.0)
		fstat=0
		mf=0.0
		DO q2=1, max_cell_read

			if(q2>max_cell_read) exit

			r(:)=x(:)+cube(q2,:)
			rr=dble(cube(q2,:))
			ra=(rr(1))**2.0 + (rr(2))**2.0 + (rr(3))**2.0
			ra=sqrt(ra)+1.!*sqrt(3.0)

			if(ra<radhii) then
			xhii=1.0
			else
			dat_overlap(5,q)=real(q2)
			exit
			end if

			call check_boundary(r, n_cell, r1, r2, r3)

			arr(r1,r2,r3)=arr(r1,r2,r3) + xhii
		

		END DO
	END DO
end if



!!round two partially ionised region..
if((myoffset >0) .and. (myonset>=myoffset)) then

	DO q=myoffset,myonset
		x(1:3)=int(dat_overlap(1:3,q)) 
		radhii=(3.0/4.0/pi*dat_overlap(4,q))**(1.0/3.0)
		den=parameter1d_array(2,q)
		fstat=0
		mf=0.0

		mfexp=(dat_overlap(4,q)-dat_overlap(5,q)+1.d0)

!DO
		DO q2=int(dat_overlap(5,q)), max_cell_read

			if(q2>max_cell_read) then
			write(*,*) 'Overshooting array..!!'
			 exit
			end if

			r(:)=x(:)+cube(q2,:)
			rr=dble(cube(q2,:))
			ra1=(rr(1))**2.0 + (rr(2))**2.0 + (rr(3))**2.0
			ra=sqrt(ra1)+1.!*sqrt(3.0)

			if(ra<radhii) then
			xhii=1.0
			else
			if((ra-radhii)<1.0*sqrt(3.0)) then

			!ra=sqrt(ra1)+1.!*sqrt(3.0)
			xhii=(radhii-dble(int(radhii)))/(ra-dble(int(radhii)))

			else
			xhii=0.0
			exit
			end if
			end if

			if(xhii<0.) xhii=0.0
			if(xhii>1.0) xhii=1.0
			call check_boundary(r, n_cell, r1, r2, r3)

			arr(r1,r2,r3)=arr(r1,r2,r3) +real(xhii/(dble(over_density(r1,r2,r3)))/den)
			mf=mf+xhii!/(dble(over_density(r1,r2,r3)))/den	
			if(mf>mfexp)	then
			fstat=1
			exit
			end if			

		END DO



	END DO
end if



	mysum=mf/(dble(nn))**3.0 

	
	!write(*,*) 'exit visual1, rank, contribution mf', myid, mysum

!if(myid==1) write(*,*) 'max min xhi: visual4', maxval(arr), minval(arr)


END SUBROUTINE visual_overlap4


subroutine xhii_xhi_binary(arr, arrhibi)
	use param
	implicit none
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr, arrhibi
	integer :: i, j, k, nn

	nn=n_cell



	DO i=1, nn
	DO j=1, nn
	DO k=1, nn
	if(arr(i,j,k)>1.0) arr(i,j,k)=1.0
		if(arr(i,j,k)==1.0) then
			arr(i,j,k)=0.0
			arrhibi(i,j,k)=0.0
		else
			arrhibi(i,j,k)=1.0
			arr(i,j,k)= 1.0-arr(i,j,k)

		end if
	END DO
	END DO
	END DO


end subroutine xhii_xhi_binary


!! This one will generate the HII bubbles with overlap. photon not conserved.
!! This also creat partial ionization fraction of the HII bubble is very small.
!!Output array nhiiarr: xhii with overlap.. Maxval(nhiiarr)--> max number of overlaps, may be larger than 1.
!!
SUBROUTINE visual_overlap1(myoffset, myonset, myid, mysum, nhiiarr)
	use param
	implicit none
	integer, intent(in) :: myoffset, myonset, myid
	real(8), intent(out) :: mysum
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: nhiiarr
	INTEGER :: q, q2, r1, r2, r3, nn, cells, con
	INTEGER,DIMENSION(3):: x, r
	real(8) ::  mf, extra

	con=0
	mf=0.0
	extra=0.0
	nn=n_cell
	nhiiarr=0.0

	if((myoffset >0) .and. (myonset>=myoffset)) then
		DO q=myoffset,myonset
			x(1:3)=int(dat_overlap(1:3,q)) 
			cells=int(dat_overlap(4,q))+1
			con=con + cells

			DO q2=1,cells
				r(:)=x(:)+cube(q2,:)
				call check_boundary(r, n_cell, r1, r2, r3)
				if(q2<cells) then
					nhiiarr(r1,r2,r3)=nhiiarr(r1,r2,r3) + 1.0
					mf=mf+over_density(r1,r2,r3)		
				else
					extra= dat_overlap(4,q)-dble(cells-1)	!!This the the partially ionization fraction
					nhiiarr(r1,r2,r3)=nhiiarr(r1,r2,r3)+real(extra)
					mf=mf+extra*over_density(r1,r2,r3)
				end if						
			END DO
		END DO
	end if

	mysum=mf/(dble(nn))**3.0 
END SUBROUTINE visual_overlap1




!!This calculate the unused photons for each halo.
!! Maxval(nhiiarr)--> max number of overlaps, may be larger than 1.
!	
SUBROUTINE visual_overlap2(myoffset, myonset, myid, nhiiarr)
	use param
	IMPLICIT NONE
	integer, intent(in) :: myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(in) :: nhiiarr
	INTEGER:: q, q2, r1, r2, r3, nn, cells
	INTEGER,DIMENSION(3):: x, r
	real(8)::loss, over
	real(4)::  rho, nhii, rad

	over=0.0
	nn=n_cell
	if((myoffset >0) .and. (myonset>=myoffset)) then
		DO q=myoffset,myonset
			x(1:3)=int(dat_overlap(1:3,q)) 
			cells=int(dat_overlap(4,q))+1
			loss=0.0

			DO q2=1,cells
				r(:)=x(:)+cube(q2,:)
				call check_boundary(r, n_cell, r1, r2, r3)
				nhii=nhiiarr(r1,r2,r3)
				rho=over_density(r1,r2,r3)
				if(nhii> 1.0) then
					loss=loss + dble((1.0-1.0/nhii) *rho) !!account for the loss photons by mass fraction
					over=over + dble((1.0/nhii) *rho) !! contribution in overlap part
				end if
			END DO
			dat_overlap(5,q) = real(loss)
		END DO
		dat_overlap(6,myoffset:myonset)=dat_overlap(4,myoffset:myonset)
	end if

END SUBROUTINE visual_overlap2

!!This one redistribute the unused photons.
!! arr input as xhii, output as xhi 
SUBROUTINE visual_overlap3(myoffset, myonset, myid, arr)
	use param
	IMPLICIT NONE
	INTEGER, intent(in) :: myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr
	INTEGER::q,r1,r2,r3,nn,cells,i,j,k,fstat, imfp, nnhalo, p
	INTEGER,DIMENSION(3)::x, r
	real(8)::extra,tot, xhii, avail, rad, rcell
	real(4), dimension(6) :: arr6
	real(4), dimension(5) :: arr5




	rcell=R_mfp*(1.d0+comm_redshift)*dble(n_cell)/box/hlittle	!!scale of mean free path in grid unit
	imfp=int(4.d0/3.d0*pi*rcell**3.d0)

	nn=n_cell

!! Now overlap is calculated, so make the dummy_nharr arr such that reionized reion is <=1, otherwise 0

	DO i=1, nn
	DO j=1, nn
	DO k=1, nn
	if(arr(i,j,k)>1.0) arr(i,j,k)=1.0	!!xhii
	END DO
	END DO
	END DO

nnhalo=num_halo



	nn=n_cell
	extra=0.0

	DO 
	fstat=0
		DO q=1,nnhalo

			if(dat_overlap(5,q)>0.05) then
				x(1:3)=int(dat_overlap(1:3,q)) 
				if(int(dat_overlap(4,q))>=1) then

				DO i=int(dat_overlap(4,q)),spdim
					r(:)=x(:)+cube(i,:)
					call check_boundary(r, n_cell, r1, r2, r3)

					!!mfp
					if(i>=imfp) then
						dat_overlap(5,q)=0.0
						exit
					end if

					dat_overlap(4,q)=dat_overlap(4,q)+1.0
					if(int(dat_overlap(4,q))>=spdim) then
						dat_overlap(5,q)=0.0
						!write(*,*) 'Overshooting...arr(q),spdim',spdim
						exit
					end if

					xhii=arr(r1,r2,r3)
					if(xhii<1.0) then
						avail=(1.0-xhii)*over_density(r1,r2,r3)
						if(dat_overlap(5,q)>avail) then
							arr(r1,r2,r3)=1.0
							dat_overlap(5,q)=dat_overlap(5,q)-avail
						else
							dat_overlap(5,q)=0.0
							arr(r1,r2,r3)=avail/over_density(r1,r2,r3)+xhii
							fstat=1
						end if
						exit
					end if	

				END DO !i
				end if
else
DO p=nnhalo, q+1, -1
	if(dat_overlap(5,p)>0.05) then
	arr6=dat_overlap(:,p)
	dat_overlap(:,p)=dat_overlap(:,q)
	dat_overlap(:,q)=arr6

arr5=parameter1d_array(:,p)
parameter1d_array(:,p)=parameter1d_array(:,q)
parameter1d_array(:,q)=arr5
	nnhalo=nnhalo-1
	exit
	else
	nnhalo=nnhalo-1
	end if
END DO
if(nnhalo<1) nnhalo=1
!!
			end if
	
		END DO !q

!write(*,*) 'RAGHU', nnhalo

		if(fstat==0) then
			exit
		end if

	END DO

write(*,*) 'max min xhi: visual3', maxval(arr), minval(arr)

	!write(*,*) 'Overlap3 exit : rank :vol frac',myid,1.0-sum(arr(:,:,:))/(dble(n_cell))**3.0


END SUBROUTINE visual_overlap3



SUBROUTINE visual_overlap3_mod(myoffset, myonset, myid, arr, fstat)
	use param
	IMPLICIT NONE
	INTEGER, intent(in) :: myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr	!!xhii .. can have larger than 1
	integer, intent(out) ::fstat
	INTEGER::q,r1,r2,r3,nn,cells,i,j,k, ii, jj, cn
	INTEGER,DIMENSION(3)::x, r
	real(8)::extra,tot, xhii, avail, rad, xhii_con, xhii_collect, exs, rho

	nn=n_cell
	extra=0.0
	fstat=0
	cn=0

!!count the used photons in this step



		DO q=myoffset, myonset
			ii=int(dat_overlap(4,q))
			jj=int(dat_overlap(6,q))+1

			if(ii>=jj) then

				x(1:3)=int(dat_overlap(1:3,q)) 
				DO i=jj, ii
				r(:)=x(:)+cube(i,:)
				call check_boundary(r, n_cell, r1, r2, r3)
				xhii=arr(r1,r2,r3)
				if(xhii>2.0) then	!!for somereason the extra area is getting min value 2 instead of 1
				rho=over_density(r1,r2,r3)
				dat_overlap(5,q)=dat_overlap(5,q)+ (1.0-1.0/xhii)*rho

				end if
				END DO

			end if


		END DO



	exs=sum(dat_overlap(5,myoffset:myonset))

	if(exs>0.5d0) then

		DO q=myoffset, myonset

!!redistribute

			if(dat_overlap(5,q)>0.05) then
			if(int(dat_overlap(4,q))>=1) then
				xhii_con=dat_overlap(5,q)
				xhii_collect=0.0
				x(1:3)=int(dat_overlap(1:3,q)) 
				dat_overlap(6,q)=dat_overlap(4,q)	!!saving the last radius

				DO i=int(dat_overlap(4,q)),spdim
					r(:)=x(:)+cube(i,:)
					call check_boundary(r, n_cell, r1, r2, r3)



					dat_overlap(4,q)=dat_overlap(4,q)+1.0
					if(int(dat_overlap(4,q))>=spdim) then
						dat_overlap(5,q)=0.0
						dat_overlap(6,q)=dat_overlap(4,q)
						exit
					end if

					xhii=arr(r1,r2,r3)
					if(xhii<1.0) then
						avail=(1.0-xhii)*over_density(r1,r2,r3)
						if(dat_overlap(5,q)>avail) then
							arr(r1,r2,r3)=1.0
							dat_overlap(5,q)=dat_overlap(5,q)-avail
							xhii_collect=xhii_collect+avail
						else
							dat_overlap(5,q)=0.0
							arr(r1,r2,r3)=avail/over_density(r1,r2,r3)+xhii
							fstat=1
							xhii_collect=xhii_con+1.0
						end if
						if(xhii_collect>xhii_con) exit

					end if	

				END DO !i
		end if
		end if
	
		END DO !q

	else
		dat_overlap(6,myoffset:myonset)=dat_overlap(4,myoffset:myonset)
	END if


END SUBROUTINE visual_overlap3_mod

subroutine det_mod_mass_fx(m)
use param
implicit none
real(8), intent(inout) :: m
integer :: i
real(8) :: mm

m=m*param_fx/param_fx_fix_1d
if(m < mass1d_min) then
	m=mass1d_min
elseif(m > mass1d_max) then
	m=mass1d_max
else
	DO i=1, n_ifront
	mm=ifront_data(i,5)
	if(mm>=m) then
	m=mm
	exit
	end if
	END DO
end if

end subroutine det_mod_mass_fx

!!This subroutine include the partial ionization regions using the RT profiles
!! arr return is xhii
subroutine visual_nh_tk(myoffset, myonset, myid, arr)
	use param
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr

	integer::n1, n2, num, ii, jj, i, j, k
	integer::i1, i2, fstat, nn, r1, r2, r3, q,  npp, inc
	integer,dimension(3)::x, r

	real(8)::T0, Tg, volc, sim_unit, rad, m, z, d, t, xhii, bb, xhii_box, sim_unit1, xhii_boxp, cnt, expt_hii, tot_hii, rad_hii1d, xhii_1d, incv, rad_1d, rad_mod, rad_hii, fx, al
	character(len=180) ifile
  	character(7) :: z_s

	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	real(4), dimension(:,:), allocatable::arr_bi
	real(4), dimension(:), allocatable::read_dat, av_read_dat

arr=0.0

if((myoffset >0) .and. (myonset>=myoffset)) then

	nn=n_cell
	sim_unit=box/dble(nn)/(1.0+comm_redshift)/hlittle !mpc physical


	DO i2=myoffset,myonset !halo
		x(:)=int(dat_overlap(1:3,i2)) 
		

		m=parameter1d_array(1,i2)
call det_mod_mass_fx(m)
		d=parameter1d_array(2,i2)
		z=parameter1d_array(3,i2)
		t=parameter1d_array(4,i2)
		fx=param_fx_fix_1d
		al=param_al_fix_1d


		n1=int(dat_overlap(4,i2))+1
		rad_Hii=(dble(n1)/1.333/pi)**(1./3.)	!!This is the radius of the overllaped big bubble

!write(*,*) rad_Hii


		call input_rt_file(m, d, z, t, source_sed, fx, al, ifile)

!ifile='/disk/dawn-1/gragh/grizzly.k/grizzly_rt/1dprof_sed3/m0.100E+09_d1.000_sed3_fx0.050_al1.500_z20.000_t10.bin'
		!open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='binary')
		open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='unformatted')

		if (fstat /= 0) then
			write(*,*) 'error opening file',ifile
		endif
		read(14) ii,jj
		allocate(arr_bi(ii,jj), read_dat(jj), av_read_dat(jj))
		read(14) arr_bi
		close(14)


!!determine the radius at which the xhii_cut 
		sim_unit1=sim_unit*(1.0+comm_redshift)/(1.0+z)
		DO i=1,ii
			read_dat=arr_bi(i, :)
!write(*,*) read_dat
			xhii=1.0-read_dat(2)
			if(xhii < (1.0-xhii_cut)) then
				rad_hii1d=arr_bi(i,1)/sim_unit1*1d-3	!!This is the radius after which the 1d profile have to put
				volc=1.333*pi*(rad_hii1d**3.0)
				npp=int(volc)	!!This is 

!write(*,*) read_dat, npp
				exit
			end if
		END DO

!stop

		read_dat=0.0


		av_read_dat=0.0
		cnt=0.0
		fstat=0

		DO i1=1,ii-1   !grid

			read_dat=arr_bi(i1, :)

			sim_unit1=sim_unit*(1.0+comm_redshift)/(1.0+z)	!!correction for the 1d z grid
			rad_1d=arr_bi(i1,1)/sim_unit1*1d-3	!!This is the radius frm the 1d profile

Rad_mod=rad_1d+rad_Hii-rad_hii1d	!!this is the distance after adding the overlapped redius

			xhii_1d=1.0-read_dat(2)
			if(xhii_1d < 1d-5) exit  !!1e-6 for fine distribution of nh,tk !1e-5 for speed

			if(xhii_1d < (1.0-xhii_cut)) then
				av_read_dat= av_read_dat + read_dat
				cnt=cnt + 1.0

			end if

			volc=1.333*pi*(rad_mod**3.0)
			n2=int(volc)

!rad=rad + rad_hii !!this is the radius beyonf the overlapped hii bubble

			tot_hii=0.0
			if(n2>=n1 .and. xhii_1d < (1.0-xhii_cut)) then	

!!write(*,*) n1, n2, read_dat
!stop

				xhii_1d=1.0-av_read_dat(2)/cnt
				expt_hii=dble(n2-n1+1)*xhii_1d*d!*100.	!!summed mass fraction!! This will dilute with radius.. as there will be bigger HII region due to overlap

					DO q=n1,max_cell_read
						!rad=(dble(q)/1.333/pi)**(1./3.)	!!this is the radius at which we are considering partial ionization
						xhii=xhii_1d*(rad_1d/rad_mod)**2.	!!This is the diluted xhii due to far away distance

!write(*,*) xhii, xhii_1d, rad_hii,rad_hii1d,rad_1d,rad_mod
!stop
						if(q>max_cell_read) exit

						r(:)=x(:)+cube(q,:)
						call check_boundary(r, n_cell, r1, r2, r3)
						xhii_box=arr(r1,r2,r3)


						xhii_boxp=xhii_box+ xhii*(1.0-xhii_box)*d/dble(over_density(r1,r2,r3))
						if(xhii_boxp >1.0) xhii_boxp=1.0




						if(arr(r1,r2,r3) < real(xhii_boxp)) then
 						arr(r1,r2,r3) =  real(xhii_boxp)
if(r1==1 .and. r2==1 .and. r3==1) write(*,*) xhii_boxp, xhii_box
						tot_hii=tot_hii+(xhii_boxp-xhii_box)*dble(over_density(r1,r2,r3))
						end if



						if(tot_hii>=expt_hii) then
						fstat=0
						exit
						end if
					END DO
					n1=q+1

!incv=4.d0*pi*rad_mod**2.d0/5.d0	!delr =1/2 grid
!if(incv>1.d0) then
!inc=int(incv)
!else
!inc=1
!end if
!				npp=n2+inc

				av_read_dat=0.0
				cnt=0.0
			end if

		END DO

	
		deallocate(arr_bi, read_dat, av_read_dat)		

!stop

	END DO !halo



end if

end subroutine visual_nh_tk


!!This subroutine include the partial ionization regions using the RT profiles



!!This subroutine include the partial ionization regions using the RT profiles
!! arr return is xhii
subroutine visual_nh_tk_incorrect(myoffset, myonset, myid, arr)
	use param
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::myoffset, myonset, myid
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arr

	integer::n1, n2, num, ii, jj, i, j, k
	integer::i1, i2, fstat, nn, r1, r2, r3, q, rad_hii, npp
	integer,dimension(3)::x, r

	real(8)::T0, Tg, volc, sim_unit, rad, m, z, d, t, xhii, bb, xhii_box, sim_unit1, xhii_boxp, cnt, expt_hii, tot_hii, rad_hii1d, xhii_1d
	character(len=180) ifile
  	character(7) :: z_s

	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	real(4), dimension(:,:), allocatable::arr_bi
	real(4), dimension(:), allocatable::read_dat, av_read_dat

arr=0.0

if((myoffset >0) .and. (myonset>=myoffset)) then

	nn=n_cell
	sim_unit=box/dble(nn)/(1.0+comm_redshift)/hlittle !mpc physical


	DO i2=myoffset,myonset !halo
		x(:)=int(dat_overlap(1:3,i2)) 
		

		m=parameter1d_array(1,i2)
		d=parameter1d_array(2,i2)
		z=parameter1d_array(3,i2)
		t=parameter1d_array(4,i2)


		n1=int(dat_overlap(4,i2))+1
		rad_Hii=(dble(n1)/1.333/pi)**(1./3.)	!!This is the radius of the overllaped big bubble


		call input_rt_file(m, d, z, t, source_sed, param_fx, param_al, ifile)
		!open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='binary')
		open(unit=14,file=trim(ifile), status='old',iostat=fstat,form='unformatted')

		if (fstat /= 0) then
			write(*,*) 'error opening file',ifile
		endif
		read(14) ii,jj
		allocate(arr_bi(ii,jj), read_dat(jj), av_read_dat(jj))
		read(14) arr_bi
		close(14)

!!determine the radius at which the xhii_cut 
		sim_unit1=sim_unit*(1.0+comm_redshift)/(1.0+z)
		DO i=1,ii
			read_dat=arr_bi(i, :)
			xhii=1.0-read_dat(2)
			if(xhii < (1.0-xhii_cut)) then
				rad_hii1d=arr_bi(i,1)/sim_unit1*1d-3	!!This is the radius after which the 1d profile have to put
				volc=1.333*pi*(rad_hii1d**3.0)
				npp=int(volc)	!!This is 
				exit
			end if
		END DO



		read_dat=0.0


		av_read_dat=0.0
		cnt=0.0
		fstat=0

		DO i1=1,ii-1   !grid

			read_dat=arr_bi(i1, :)

			sim_unit1=sim_unit*(1.0+comm_redshift)/(1.0+z)	!!correction for the 1d z grid
			rad_hii1d=arr_bi(i1,1)/sim_unit1*1d-3	!!This is the radius frm the 1d profile

			xhii_1d=1.0-read_dat(2)
			if(xhii_1d < 1d-6) exit  !!1e-6 for fine distribution of nh,tk !1e-5 for speed

			if(xhii_1d < (1.0-xhii_cut)) then
				av_read_dat= av_read_dat + read_dat
				cnt=cnt + 1.0

			end if

			volc=1.333*pi*(rad_hii1d**3.0)
			n2=int(volc)

!rad=rad + rad_hii !!this is the radius beyonf the overlapped hii bubble

			tot_hii=0.0
			if(n2>=npp .and. xhii_1d < (1.0-xhii_cut)) then	

				xhii_1d=1.0-av_read_dat(2)/cnt
				expt_hii=dble(n2-npp+1)*xhii_1d*d	!!summed mass fraction!! This will dilute with radius.. as there will be bigger HII region due to overlap

					DO q=n1,max_cell_read
						rad=(dble(q)/1.333/pi)**(1./3.)	!!this is the radius at which we are considering partial ionization
						xhii=xhii_1d*(rad_hii1d/rad)**2.	!!This is the diluted xhii due to far away distance
						if(q>max_cell_read) exit

						r(:)=x(:)+cube(q,:)
						call check_boundary(r, n_cell, r1, r2, r3)
						xhii_box=arr(r1,r2,r3)


						xhii_boxp=xhii_box+ xhii*(1.0-xhii_box)*d/dble(over_density(r1,r2,r3))
						if(xhii_boxp >1.0) xhii_boxp=1.0
						tot_hii=tot_hii+(xhii_boxp-xhii_box)*dble(over_density(r1,r2,r3))



						if(arr(r1,r2,r3) < real(xhii_boxp)) arr(r1,r2,r3) = real(xhii_boxp)

						if(tot_hii>=expt_hii) then
						fstat=0
						 exit
						end if
					END DO
					n1=q+1
				npp=n2+1

				av_read_dat=0.0
				cnt=0.0
			end if

		END DO

	
		deallocate(arr_bi, read_dat, av_read_dat)		


	END DO !halo



end if

end subroutine visual_nh_tk_incorrect

!!This provides the effective temperature maps.. this can be different than the actual temperature maps.. due to resolution issues.
!! 
subroutine gen_tkmap(myid, z, arrxhi, arrtk)
	use param
	implicit none
	integer, intent(in) :: myid
	real(8), intent(in) :: z
	real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arrxhi
	real(4), dimension(n_cell, n_cell, n_cell), intent(out) :: arrtk
	integer :: i, j, k, ii
	real(8) :: xhii

	call read_reion(myid, z)	!!this will initialize corel array


	arrtk=(tcmb*(1.0+z)*(1.0+z)/(z_dec))

	corel(1,1)=1.0
	DO k=1,n_cell
	DO j=1,n_cell
	DO i=1,n_cell

		if(arrxhi(i,j,k) <0.0) then
		write(*,*)'Error'
		 stop
		end if

		if(arrxhi(i,j,k) >1.0) then
		arrxhi(i,j,k)=1.0
		end if

		xhii=1.0-arrxhi(i,j,k)

		DO ii=1,n_corel
			if(corel(ii,1) <= xhii) then
				arrtk(i,j,k)=corel(ii,2)
				exit
			end if
		END DO

	END DO
	END DO
	END DO


	write(*,*) 'max min xhi', maxval(arrxhi), minval(arrxhi)
	write(*,*) 'max min tk', maxval(arrtk), minval(arrtk)


end subroutine gen_tkmap

subroutine bound_checkxhi(arrxhi)
use param
implicit none
real(4), dimension(n_cell, n_cell, n_cell), intent(inout) :: arrxhi
integer :: i, j, k


DO i=1, n_cell
DO j=1, n_cell
DO k=1, n_cell
if(arrxhi(i,j,k)>1.0) arrxhi(i,j,k)=1.0
if(arrxhi(i,j,k)<0.0) arrxhi(i,j,k)=0.0
END DO
END DO
END DO

write(*,*) 'Bound checked', maxval(arrxhi), minval(arrxhi)
end subroutine bound_checkxhi



