MODULE param

!!1D CODE PARAMETERS
	real(8), parameter:: hplanck_erg = 6.6260695729E-27
	real(8), parameter:: hplanck_ev = 4.13566751691E-15
	real(8), parameter:: erg_ev = 6.2415E+11
	real(8), parameter:: E_HI = 13.60	
	real(8), parameter:: E_lyalpha = 10.20		
	real(8), parameter:: E_HeI = 24.590		
	real(8), parameter:: E_HeII = 54.420	
	real(8), parameter:: E_max = 10000.00	
	real(8), parameter:: Megayear = 3.1536E+13
	real(8), parameter:: Megaparsec = 3.08568025E+24                 !megaparsec in cm
	real(8), parameter:: kiloparsec = 3.08568025E+21			!kiloparsec in cm
	real(8), parameter:: angstrom_cm = 1.d-8
	real(8), parameter:: pi = 3.14159265359
	real(8), parameter:: C_light = 3.0e+10 
	real(8), parameter:: lambda_alpha = 1216.0	!!In A
	real(8), parameter:: lambda_beta = 1026.0	!!In A
	real(8), parameter:: lambda_hi = 912.0	!!In A
  	real(8), parameter :: grav_const = 6.672d-8 !cm^3/g/s^2
  	real(8), parameter :: Gkgs = grav_const/1.0d3
  	real(8), parameter :: M_sun = 1.989d33
	real(8), parameter :: m_proton = 1.6726d-24	!!gram
	real(8),parameter:: Kboltz_ev=8.62911325d-5 
	real(8),parameter:: L_solar_erg=3.839d+33 !erg/s		

!! COSMOLOGICAL PARAMETERS .. suman da paper : WMAP5

	real(8), parameter :: omega_l = 0.73
	real(8), parameter :: omega_m = 1.0 - omega_l
	real(8), parameter :: omega_k = 1.0 - omega_l - omega_m
	real(8), parameter :: omega_r = 5.0E-5
	real(8), parameter :: omega_b = 0.044	!!not 0.046
	real(8), parameter :: omega_ch = 0.7
	real(8), parameter :: bias = 1
	real(8), parameter :: power_index = 2.0 
	real(8), parameter :: hlittle = 0.7
	real(8), parameter :: Ho = (hlittle*3.2407e-18) ! s^-1 at z=0 
	real(8), parameter :: hydrogen_frac = 0.752
	real(8), parameter :: f_he = 1.0 - hydrogen_frac
	real(8), parameter :: z_dec = 151.0
	real(8), parameter :: tcmb = 2.73
	real(8), parameter :: h0kgs = hlittle*3.2407e-18 
	real(8), parameter :: rho_ckgs = (3.0*h0kgs**2)/(8.0*pi*Gkgs)

!!!SIMULATION PARAMETERS


	real(8), parameter :: box = 114.0 !!box length cMpc/h unit
	integer, parameter :: n_cell =  256	!!Number of grid in the simulation box 
	integer, parameter :: nc = 6144	!!This helps to read calculate density !!only for C2RAY

!!PARAMETERS

	real(8) :: param_zeta, param_Mmin, param_fx, param_al


	integer, parameter :: model_sc = 1	!!1 ts>>tcmb, 2 ts=tk, 3 ts self
	integer, parameter :: source_sed = 6

	real(8), parameter :: fix_age = 100.0
	real(8), parameter :: param_zeta_fix = 1.d0!5.d0*1.479d-3! 1.d0
	real(8), parameter :: param_Mmin_fix = 2.2d9
	real(8), parameter :: param_fx_fix_1d = 0.0d0
	real(8), parameter :: param_al_fix_1d = 0.0d0

	integer, parameter :: nparam1d = 7

	logical, parameter :: manual_terminate = .true.
	integer, parameter :: nterminate = 3

!!RELATED TO SIMULATION

	integer, parameter :: parameter_space = 0	!!set 1 if want to explore parameter space, else set 0
	integer, parameter :: track_age = 1 	!!set 1 if want to track the age of halo, otherwise will set age 10 Myr
	integer, parameter :: genxhimap = 1	!!set 1 if want to generate xhi maps, 0 otherwise
	integer, parameter :: useproilexhi = 0	!!set 1 if want to use Tb profiles to generate xhi maps, 0 otherwise
	integer, parameter :: gentkmap = 0	!!set 1 if want to generate TK maps, 0 otherwise
	integer, parameter :: genalphamap = 0	!! set 1 if want to generate alpha maps, 0 otherwise
	integer, parameter :: gentbmap = 0	!!set 1 if want to generate Tb maps, 0 otherwise
	integer, parameter :: useanalytic = 0	!!set 1 if one to use stromgen analititcal formula, else put 0
	integer, parameter :: genps = 0		!!set 1 if want to generate power spectrum, else put 0.. for tb
	integer, parameter :: genpsxx = 0		!!set 1 if want to generate power spectrum, else put 0.. for xx
	integer, parameter :: genpsxd = 0		!!set 1 if want to generate power spectrum, else put 0.. for xd
	integer, parameter :: selected_reion = 1	!!if 1 then only consider selected reionization redshifts, 0 will create entire history

	integer, parameter :: writexhimap = 1
	integer, parameter :: writetkmap = 0	
	integer, parameter :: writealphamap = 0	
	integer, parameter :: writetbmap = 0

	logical,parameter:: visual_multicore = .false.	!!this will apply parallel verlap distribution

	logical,parameter::restart_para_est=.false.	!!Set true if want to restart the simulation from certain redshift 
	INTEGER,PARAMETER::restart_para_est_id=360	!!this is the number of jobs which are completed	

	integer, parameter :: gr_members_max = 5	!!number of cores per task


!!ARRAYS used in the simulation..

	integer, parameter :: max_input=300	!!Number of redshift bins
	real(8), dimension(max_input) :: z_checkpoint, age_checkpoint, z_selectpoint


	integer::max_halos
	real(4), dimension(:, :), allocatable:: dat_overlap, parameter1d_array	!!This carry all the information about the haloes


	integer, parameter::n_corel = 30! int(c_light*50.0*megayear/kiloparsec/20.0)	!!This is required to carry the Tk, xHI correlation
 	real(8), dimension(n_corel, 2):: corel


	integer, parameter :: ndum = int(150.d0/box*dble(n_cell))	!!max escpae fraction of the array=150/h mpc
	integer, parameter :: ncube=min(n_cell, ndum)
	integer, parameter:: max_cell_read=(ncube-1)**3	!!This is the dimension of array cube, which helps to create the sphere around the sources.
	integer, dimension(max_cell_read, 3)::cube !!This array will be used for creating spherical bubbles


	integer:: n_ifront, check_min
	real(4), dimension(:, :), allocatable:: ifront_data

!! SIMULATION PARAMETERS

 	character(*), parameter :: comm_path    = '/disk/dawn-1/gragh/nbody_data/nb_l114_ncell256_nc_6144/' 
 	character(*), parameter :: denisty_path    = comm_path//'density_fields/nc256_halos_removed/'
 	character(*), parameter :: halolist_path = comm_path//'sourcelist/source_2.2e9/'
 	character(*), parameter :: output_path    = '../output_grizzly/' 
	character(*), parameter :: output_1d      = '/disk/dawn-1/gragh/grizzly.k/grizzly_rt/1dprof_sed6/'!'/disk/dawn-1/gragh/grizzly.k/grizzly_rt/1dprof_sed3/'
 	character(*), parameter :: checkpoints ='/disk/dawn-1/gragh/grizzly.k/grizzly_sim/inputs/reion' 
 	character(*), parameter :: selectpoints ='/disk/dawn-1/gragh/grizzly.k/grizzly_sim/inputs/select_z' 
 	character(*), parameter :: cellpath='/disk/dawn-1/gragh/grizzly.k/grizzly_sim/inputs/'
 	character(*), parameter :: pspath = '/disk/dawn-1/gragh/grizzly.k/grizzly_sim/inputs/'	!!put the SED files here
 

!!Information for selecting the 1D prfile


	real(8) :: mass1d_max, mass1d_min, delta1d_max, delta1d_min, age1d_max, age1d_min, z1d_max, z1d_min, fx1d_min, fx1d_max, al1d_min, al1d_max

	real(8), parameter :: xhii_cut = 0.5
	character(180) :: filename_comm


!! Large Arrays used in the simulation 

	integer, dimension(n_cell, n_cell, n_cell):: haloagenew
	real(4), dimension(n_cell, n_cell, n_cell):: over_density, matrix_nhI, matrix_tk, matrix_alpha, halomassnew


	real(8):: Mhalo_tot, comm_redshift, gap_timestep, gap_1sttimestep, Total_sed_energy, Total_sed_photon, mhalo_tot_L1, mhalo_tot_L2, mhalo_tot_L, global_expt_mf
	integer :: spdim, num_halo, rank, numtasks, num_checkpoints, redshift_loop_index, num_selectpoints


	integer, parameter :: Nbins_ps = 20 !Division of the k bin



!! Large Arrays used in the simulation 

	real(8), parameter :: g_gamma_c2ray_L = 25.0
	real(8), parameter :: g_gamma_c2ray_S = 0.0
	real(8), parameter :: f_esc_hi = 0.1
	real(8), parameter :: f_esc_lya = 1.0	!!This should be 1.. 0.1 to incorporate f_star!!
	real(8), parameter :: f_lya = 3.0
	real(8), parameter :: n_gamma_hi = 5.32d41	!!number of ionizing photons per sec per solar mass halo
	logical, parameter :: stellar_sed = .true.

	integer, parameter :: nsed = 141	!50
	real(8), dimension(nsed) :: lambdaarr, lumarr

	real(8):: time_since, lya_redshift
	integer :: r_grid_lalpha_min, r_grid_lalpha_max, imax_comm, imin_comm


	character(180), parameter :: halofilename = '-coarsest_sources.dat'
	character(180), parameter :: densityfilename = 'n_all.dat'
	character(180), parameter :: velfilename = ''
	character(180), parameter :: xhifilename = 'xhi.bin'
	character(180), parameter :: tkfilename = 'tk.bin'
	character(180), parameter :: alphafilename = 'lyalpha.bin'
	character(180), parameter :: tbfilename = 'tb.bin'
	character(180), parameter :: pskfilename = 'PSK.dat' 
	character(180), parameter :: pskxxfilename = 'PSKxx.dat' 
	character(180), parameter :: pskxdfilename = 'PSKxd.dat' 

	integer, parameter :: n_param_sim = 4
	integer, parameter :: nn_zeta=1
	integer, parameter :: nn_Mmin=1
	integer, parameter :: nn_fx=4
	integer, parameter :: nn_al=1

	integer, parameter :: n_jobs = nn_zeta*nn_mmin*nn_fx*nn_al
	real(8), dimension(n_jobs,n_param_sim) :: job_arr

	real(8), dimension(nn_zeta)::zeta_arr
	real(8), dimension(nn_mmin)::mmin_arr
	real(8), dimension(nn_fx)::fx_arr
	real(8), dimension(nn_al)::al_arr

        real(8), dimension(2), parameter::zeta_points=(/1.5d0,1.5d0/)
        real(8), dimension(2), parameter ::mmin_points=(/9.34d0,9.34d0/) !!10^ !!from center of the box in arcmin
        real(8), dimension(2), parameter::fx_points=(/0.d0,2.d0/)
        real(8), dimension(2), parameter::al_points=(/1.5d0,1.5d0/)

real(8), parameter :: R_mfp = 10.d0	!physical comoving mean free path


END MODULE param



