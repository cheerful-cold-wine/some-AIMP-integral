!--------------------------------------------------------------------
! Ab Initio Model Potential
! M2 Coulmb: cGTO, s-type gaussian
! integral
!--------------------------------------------------------------------
subroutine aimp_M2int(nftmax,MxAng,ndimlx,ndimrx,               &
                      Kmin,Kmax,la,lb,nSubAB,nprima,nprimb,     &
                      nx,ny,nz,                                 &
                      CoordA,CoordB,ExpCoefA,ExpCoefB,          &
                      maxnumM2,nM2,CoordM2,M2exp,M2coef,gprim)
  implicit none
  integer,intent(in)  :: nftmax,MxAng,ndimlx,ndimrx
  integer,intent(in)  :: Kmin(MxAng),Kmax(MxAng)
  integer,intent(in)  :: la,lb,nSubAB,nprima,nprimb
  integer,intent(in)  :: nx(nftmax)
  integer,intent(in)  :: ny(nftmax)
  integer,intent(in)  :: nz(nftmax)
  real*8, intent(in)  :: CoordA(3),CoordB(3)
  real*8, intent(in)  :: ExpCoefA(nprima),ExpCoefB(nprimb)
  integer,intent(in)  :: maxnumM2,nM2
  real*8, intent(in)  :: CoordM2(3)
  real*8, intent(in)  :: M2exp(maxnumM2),M2coef(maxnumM2)
  !------------------------------------------------------------------
  ! Only single gx/gy/gz is useless.
  ! M2 bring many different gx/gy/gzs, sum of them are useful.
  real*8 :: gx(0:ndimlx-1,0:ndimrx-1)
  real*8 :: gy(0:ndimlx-1,0:ndimrx-1)
  real*8 :: gz(0:ndimlx-1,0:ndimrx-1)
  real*8, intent(out) :: gprim(nSubAB*nprima*nprimb)
  !------------------------------------------------------------------
  integer :: nc0,nc1,iprim,jprim,iM2,n,i,jend,j,jx,jy,jz,ix,iy,iz
  real*8  :: ABx,ABy,ABz,rAB,Cx,Cy,Cz,M2Cx,M2Cy,M2Cz,rM2C
  real*8  :: Dx,Dy,Dz,DAx,DAy,DAz
  real*8  :: pi,pi32,E_AB,E_CM,tmpn
  real*8  :: alpha1,alpha2,alpha12,alpha3,alpha123

  pi=4.d0*datan(1.d0)
  pi32=dsqrt((pi**3))

  ! (A,alpha1|M2,alpha3|B,alpha2)
  ABx=CoordA(1)-CoordB(1)
  ABy=CoordA(2)-CoordB(2)
  ABz=CoordA(3)-CoordB(3)
  rAB=ABx**2+ABy**2+ABz**2

  gprim=0.d0
  nc0=1
  nc1=1
  do jprim=1, nprimb
    alpha2=ExpCoefB(jprim)
    do iprim=1, nprima
      alpha1=ExpCoefA(iprim)
      alpha12=alpha1+alpha2

      ! Gaussian Product Theorem
      ! A,B ---1stGPT---> C
      Cx=(alpha1*CoordA(1)+alpha2*CoordB(1))/alpha12
      Cy=(alpha1*CoordA(2)+alpha2*CoordB(2))/alpha12
      Cz=(alpha1*CoordA(3)+alpha2*CoordB(3))/alpha12
      M2Cx=CoordM2(1)-Cx
      M2Cy=CoordM2(2)-Cy
      M2Cz=CoordM2(3)-Cz
      rM2C=M2Cx**2+M2Cy**2+M2Cz**2
      E_AB=dexp(-alpha1*alpha2*rAB/alpha12)

      do iM2=1, nM2
        alpha3=M2exp(iM2)
        alpha123=alpha12+alpha3
        ! C,M2 ---2ndGPT---> D
        Dx=(alpha12*Cx+alpha3*CoordM2(1))/alpha123
        Dy=(alpha12*Cy+alpha3*CoordM2(2))/alpha123
        Dz=(alpha12*Cz+alpha3*CoordM2(3))/alpha123
        DAx=Dx-Ax
        DAy=Dy-Ay
        DAz=Dz-Az
        E_CM=dexp(-alpha12*alpha3*rM2C/alpha123)

        ! Vertical/Horizontal Recursion Relation
        ! VRR
        gx(0,0)=1.d0
        gy(0,0)=1.d0
        gz(0,0)=M2coef(iM2)*E_AB*E_CM*pi32/(alpha123*dsqrt(alpha123))

        if(ndimlx-2.gt.0) then
          gx(1,0)=DAx
          gy(1,0)=DAy
          gz(1,0)=DAz*gz(0,0)
        endif
        do n=1, ndimlx-2
          tmpn=0.5d0*dble(n)/alpha123
          gx(n+1,0)=tmpn*gx(n-1,0)+DAx*gx(n,0)
          gy(n+1,0)=tmpn*gy(n-1,0)+DAy*gy(n,0)
          gz(n+1,0)=tmpn*gz(n-1,0)+DAz*gy(n,0)
        enddo
        ! HRR
        do i=1, ndimlx-1
          jend=min0(ndimrx-1,i)
          do j=1, jend
            gx(i-j,j)=gx(i-j+1,j-1)+ABx*gx(i-j,j-1)
            gy(i-j,j)=gy(i-j+1,j-1)+ABy*gy(i-j,j-1)
            gz(i-j,j)=gz(i-j+1,j-1)+ABz*gz(i-j,j-1)
          enddo
        enddo

        nc0=nc1
        do j=kmin(lb+1),kmax(lb+1)
          jx=nx(j)
          jy=ny(j)
          jz=nz(j)
          do i=kmin(la+1),kmax(la+1)
            ix=nx(i)
            iy=ny(i)
            iz=nz(i)
            ! gprim(ncartA,ncartB,nprimA,nprimB)
            ! (ncartA*ncartB)*(nprimA*nprimB)
            gprim(nc0)=gprim(nc0)+gx(ix,jx)*gx(iy,jy)*gz(iz,jz)
            nc0=nc0+1
          enddo
        enddo
      enddo ! Cycle nM2
      nc1=nc0 ! go to next prim-block of gprim
    enddo
  enddo
  !nc1=nc1-1
  return
end subroutine aimp_M2int
