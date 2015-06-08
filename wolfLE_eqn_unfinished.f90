!Based on A. Wolf 1980

program ode
implicit none
!N = # of nonlinear equations
!NN = Total # of eqtns
integer,parameter::N=3;
integer,parameter::NN=12
integer::nstep,neq,ind

real::tol,stpsze,io,x
real::y(NN),znorm(N),gsc(N),cum(N),c(24),w(NN,9)

external FCN

!
!   Initial conditions for nonlinear system
!
y(1) = 10.0
y(2) = 1.0
y(3) = 0.0

!
!   Initial conditions for linear system (Orthonormal frame)
!
do i = N+1,NN
    y(i) = 0.0
end do
do i = 1,N
    y((N+1)*i) = 1.0
    cum(i) = 0.0
end do

!
!   Integration tolerance, # of integrations steps, time per step, and I/O rate
!
type*, 'TOL,NSTEP,STPSZE,IO?'
accept*,tol,nstep,stpsze,io

!
!   Initialization for integration
!
neq = NN
x = 0.0
ind =1

do i = 1,nstep
    xend = stpsze*float(i)
    call devek(neq,fcn,x,y,xend,tol,ind,c,neq,w,ier)
    !
    !   Construct a new orthonormal basis by Gram-Schmidt method
    !
    !   Normalize first vector
    !
    znorm(1) = 0.0
    do j = 1,N
        znorm(1) = znorm(1) + y(N*j+1)**2
    end do
    znorm(1) = sqrt(znorm(1))
    do j = 1,N
        y(N*j+1) = y(N*j+1)/znorm(1)
    end do

    !
    !   Generate the new orthonormal set of vectors.
    !
    do j = 2,N
    !
    !   Generate j-1 GSR coefficients
    !
        do k = 1,j-1
            gsc(k) = 0.0
            do l = 1,N
                gsc(k) = gsc(k) + y(N*l+j)*y(N*l+k)
            end do
            !
            !   Construct a new vector
            !
            do k = 1,N
                do l = 1,j-1
                    y(N*k+j) = y(N*k+j) - gsc(l)*y(N*k+l)
                end do
            end do
            !
            !   Calculate the vector's norm
            !

