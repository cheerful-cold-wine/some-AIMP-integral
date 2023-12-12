!--------------------------------------------------------------------
! Ab Initio Model Potential
! M1 Coulmb: cGTO/r, nuclear attraction times s-type gaussian
! integral
!--------------------------------------------------------------------
subroutine aimp_M1int(nftmax,MxAng,ndimlx,ndimrx,               &
                      Kmin,Kmax,la,lb,nSubAB,nprima,nprimb,     &
                      nx,ny,nz,                                 &
                      CoordA,CoordB,ExpCoefA,ExpCoefB,          &
                      maxnumM1,nM1,CoordM1,M1exp,M1coef,gprim)
  implicit none 
  integer,intent(in)  :: nftmax,MxAng,ndimlx,ndimrx
  integer,intent(in)  :: Kmin(MxAng),Kmax(MxAng)
  integer,intent(in)  :: la,lb,nSubAB,nprima,nprimb
  integer,intent(in)  :: nx(nftmax)
  integer,intent(in)  :: ny(nftmax)
  integer,intent(in)  :: nz(nftmax)
  real*8, intent(in)  :: CoordA(3),CoordB(3)
  real*8, intent(in)  :: ExpCoefA(nprima),ExpCoefB(nprimb)
  integer,intent(in)  :: maxnumM1,nM1
  real*8, intent(in)  :: CoordM1(3)
  real*8, intent(in)  :: M1exp(maxnumM1),M1coef(maxnumM1)
  !------------------------------------------------------------------
  ! Only single gx/gy/gz is useless.
  ! M1 bring many different gx/gy/gzs, sum of them are useful.
  real*8 :: gx(0:ndimlx-1,0:ndimrx-1)
  real*8 :: gy(0:ndimlx-1,0:ndimrx-1)
  real*8 :: gz(0:ndimlx-1,0:ndimrx-1)
  real*8, intent(out) :: gprim(nSubAB*nprima*nprimb)
  !------------------------------------------------------------------
  integer, parameter :: MxRysRoot=14
  !------------------------------------------------------------------
  real*8 :: RysRt(MxRysRoot),RysWt(MxRysRoot)
  real*8 :: pi,E_AB,E_CM
  real*8 :: tmpn,fac
  real*8 :: ABx,ABy,ABz,rAB,Cx,Cy,Cz,M1Cx,M1Cy,M1Cz,rM1C,DAx,DAy,DAz
  real*8 :: alpha1,alpha2,alpha12
  integer :: nRysRoot,iRysRoot
  integer :: nc0,nc1,iprim,jprim,iM1,n,i,j,jend,ix,iy,iz,jx,jy,jz

  pi=4.d0*datan(1.d0)
  nRysRoot=(la+lb)/2+1

  ! (A,alpha1|M1,alpha3|B,alpha2)
  ABx=CoordA(1)-CoordB(1)
  ABy=CoordA(2)-CoordB(2)
  ABz=CoordA(3)-CoordB(3)
  rAB=ABx**2+ABy**2+ABz**2

  gprim=0.d0
  nc0=1
  nc1=1
  do jprim=1,nprimb
    alpha2=expcoefB(jprim)
    do iprim=1,nprima
      alpha1=expcoefA(iprim)
      alpha12=alpha1+alpha2

      ! Gaussian Product Theorem
      ! A,B ---1stGPT---> C
      Cx=(alpha1*CoordA(1)+alpha2*CoordB(1))/alpha12
      Cy=(alpha1*CoordA(2)+alpha2*CoordB(2))/alpha12
      Cz=(alpha1*CoordA(3)+alpha2*CoordB(3))/alpha12
      M1Cx=CoordM1(1)-Cx
      M1Cy=CoordM1(2)-Cy
      M1Cz=CoordM1(3)-Cz
      rM1C=M1Cx**2+M1Cy**2+M1Cz**2
      E_AB=dexp(-alpha1*alpha2*rAB/alpha12)

      do iM1=1, nM1
        alpha3=M1exp(iM1)
        alpha123=alpha12+alpha3
        ! C,M1 ---2ndGPT---> D
        E_CM=dexp(-alpha12*alpha3*rM1C/alpha123)
        ! Rys quantic
        tmpn=alpha12*alpha12*rM1C/alpha123 ! X in e^(-Xt^2)
        RysRT=0.d0
        RysWT=0.d0
        if(nRysRoot.le.3) call bdf_cvwint_rt123(tmpn,RysRT,RysWT,nRysRoot)
        if(nRysRoot.eq.4) call bdf_cvwint_root4(tmpn,RysRT,RysWT,nRysRoot)
        if(nRysRoot.eq.5) call bdf_cvwint_root5(tmpn,RysRT,RysWT,nRysRoot)
        if(nRysRoot.gt.5) call eri_RysRtsWgh(tmpn,1,RysRT,RysWT,nRysRoot)

        do iRysRoot=1,nRysRoot
          ! Vertical/Horizontal Recursion Relation
          ! VRR
          if(nRysRoot.gt.5) then
            fac  = RysRT(iRysRoot)  ! t^2 = u^2/(alpha123+u^2)
            tmpn = fac/(1.d0-fac)   ! u^2/alpha123 = t^2/(1-t^2)
            fac  = 1.d0-fac         ! 1-t^2
          else
            tmpn = RysRT(iRysRoot)  ! u^2/alpha123
            fac  = 1.d0/(1.d0+tmpn) ! alpha123/(alpha123+u^2)
          endif

          DAx=(Cx+tmpn*CoordM1(1))*fac-CoordA(1)
          DAy=(Cy+tmpn*CoordM1(2))*fac-CoordA(2)
          DAz=(Cz+tmpn*CoordM1(3))*fac-CoordA(3)
          
          gx(0,0)=1.d0
          gy(0,0)=1.d0
          gz(0,0)=2*M1coef(iM1)*E_AB*E_CM*RysWT(iRysRoot)*pi/alpha123
          if((ndimlx-2).gt.0) then
            gx(1,0)=DAx
            gy(1,0)=DAy
            gz(1,0)=DAz*gz(0,0)
          endif

          fac=alpha123*(1.d0+tmpn)
          fac=0.5d0/fac
          do n=1,ndimlx-2
            tmpn=fac*dble(n)
            gx(n+1,0)=tmpn*gx(n-1,0)+DAx*gx(n,0)
            gy(n+1,0)=tmpn*gy(n-1,0)+DAy*gy(n,0)
            gz(n+1,0)=tmpn*gz(n-1,0)+DAz*gz(n,0)
          enddo
          ! HRR
          do i=1,ndimlx-1
            jend=min0(ndimlx-1,i)
            do j=1,jend
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
            do i=kmin(la+1), kamx(la+1)
              ix=nx(i)
              iy=ny(i)
              iz=nz(i)
              ! gprim(ncartA,ncartB,nprimA,nprimB)
              ! (ncartA*ncartB)*(nprimA*nprimB)
              gprim(nc0)=gprim(nc0)+gx(ix,jx)*gz(iy,jy)*gz(iz*jz)
              nc0=nc0+1
            enddo
          enddo
        enddo ! Cycle Rys
      enddo ! Cycle nM1
      nc1=nc0 ! go to next prim-block of gprim
    enddo
  enddo
  !nc1=nc1-1
  !gprim(1:nc1)=2.d0*gprim(1:nc1)
  return
end subroutine aimp_M1int
