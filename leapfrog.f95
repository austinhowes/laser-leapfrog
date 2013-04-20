module constants

  implicit none
  complex :: pi = acos(-1.d0), i_c = cmplx(0.d0, 1.d0)

end module constants

program main
    use constants
    implicit none
    integer :: matrix_rows, matrix_cols, x, t
    double precision ss
    real :: psi_r, psi_im, dt, dx, norm
    real :: print_offset
    
    complex*16 czero,chalf,cone,ctwo,cthree,cfour 
    complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
    complex*16 a,b,c,r,gam,u,cmn,cpl,a1,b1,c1,r1,gam1,u1,psi_matrix,pot,ff,bb
    complex*16 bet,xj
    real idint, agbr, gbr, delay, del2, plato2, fwhm2, w2, e2, del1, plato1, fwhm1, w1, e1, nerg, emax, de, emin
    integer key5, key4, key3, key2, key1, ll, nn, nf,  nshape, npulse
    integer nfmax, ncmax, ntmax, nxmax, nprmax

    parameter(nxmax=381000,ntmax=100000,ncmax=61,nprmax=20000, nfmax=10)
    
    common/discr/nf,nn(nfmax),ll(nfmax)
    common/carry3/psi_matrix(0:nxmax,0:ncmax), pot(0:nxmax,0:ncmax)
      
    open(50, file='hyd7.inp', status='old')

    ! leapfrog testing file
    open(95, file="leapfrogtesting.out", status="unknown")
    open(96, file="initial.out", status="unknown")
    open(97, file="norm.out", status="unknown")
   

    ! set potential = 0 for now
    pot = 0

    !read the initial data

    read(50,*) dt, dx

    read(50,*) matrix_rows, print_offset, matrix_cols, key1, key2, key3, key4, key5

! read in quantum numbers
    if(key1 /= 0) then
        read(50,*) nf
        read(50,*) nn(1), ll(1), nn(2), ll(2), nn(3), ll(3), nn(4), ll(4), nn(5), ll(5)
!            do i = 1, nf
!                read(50,*) nn(i)
!                read(50,*) ll(i)
!                write(*,*) nn(i), ll(i)
!            end do
    end if

    if(key4 == 1) then
       read(50,*) emin, de, emax
    end if

    !nerg = idint((emax-emin)/de) + 1
    read(50,*) npulse

    if(npulse /= 0) then
        read(50,*) nshape
        read(50,*) e1,w1,fwhm1,plato1,del1

        if(npulse == 2) then
            read(50,*) e2,w2,fwhm2,plato2,del2,delay
        end if
    end if

    read(50,*) gbr, agbr

    close(50)

    call setup()
    
    do x = 0, matrix_rows
       write(96, *) x*dx, real(psi_matrix(x,0)), imag(psi_matrix(x,0)), &
            real(psi_matrix(x,0)) ** 2 + aimag(psi_matrix(x,0))**2
    end do

    !calculate second time step
      psi_matrix = psi_matrix * exp(i_c * dt / 4)
      ss = dt / dx ** 2

    print *,'starting main loop'
    
    do t = 3, ntmax
        !leapfrog propagation
        norm = 0
        do x = 1, matrix_rows - 1
           psi_r = real(psi_matrix(x,0))
           psi_im = aimag(psi_matrix(x,0))
            if (mod(t, 2) == 1) then
               ! update the real part
               psi_r = psi_r - ss * imag(psi_matrix(x + 1, 0)) &
                    - ss * imag(psi_matrix(x - 1,0)) &
                    + 2 * (ss + real(pot(x, 0))*dt) * psi_im
            else
               ! update the imaginary part
               psi_im = psi_im + ss * real(psi_matrix(x + 1,0)) &
                    + ss * real(psi_matrix(x - 1, 0)) &
                    - 2 * (ss + real(pot(x, 0))*dt) * psi_r
            end if
            psi_matrix(x,0) = cmplx(psi_r, psi_im)
            norm = norm + dx * (real(psi_matrix(x,0))**2 + aimag(psi_matrix(x,0))**2)
        end do
        ! write out norm(psi) at each time step
        write(97, *), t, norm

    end do

    ! output final wavefunction norm here
    do x = 1, matrix_rows
        write(95,*) x * dx, real(psi_matrix(x,0)), aimag(psi_matrix(x,0)), &
             real(psi_matrix(x,0))**2 + aimag(psi_matrix(x,0))**2
    end do

    contains
      subroutine setup
        ! setup psi, pot
        do x = 1, matrix_rows
           psi_matrix(x, 0) = dcmplx(2.d0 * dx * x * exp(- dx * x), 0.d0)
           pot(x, 0) = dcmplx(-1.d0 / x, 0.d0)

        end do
        psi_matrix(0,0) = cmplx(0.d0, 0.d0)
        return
      end subroutine setup
end program main
