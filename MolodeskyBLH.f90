      real*8 function MolodeskyBLH(BLH,gra,ksi,hgt,dwm,nlat,nlon,hd,dr,knd,GRS)
      !按严密积分公式计算地面及地球外部Molodensky I项
      !dr-积分半径m
!-------------------------------------------------------------
      implicit none
	integer::knd,i,j,nlat,nlon,i0,j0,ni,nj,kfn !奇异积分核函数精化参数2*kfn
	real*8::dr,gra(nlat,nlon),ksi(nlat,nlon),hgt(nlat,nlon),dwm(nlat,nlon)
	real*8::hd(6),pi,RAD,ds,mdr,tt,rr,r0,r1,r2,dwm0
	real*8::GRS(6),BLH(3),XYZ(3),rln(3),hgt0,ksi0,BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,MolodeskySgn,dl,SK,rst,gr,NFD(5)
!-----------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;kfn=4
      dwm0=CGrdPntD2(BLH(2),BLH(1),dwm,nlat,nlon,hd)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1);call BLH_XYZ(GRS,BLH,XYZ)
      mdr=rr*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/rr/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/rr/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      rst=0.d0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=hgt(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          dl=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          if(dl<mdr)then!计算奇异积分
             rst=rst+MolodeskySgn(BLH,gra,ksi,hgt,dwm,nlat,nlon,hd,knd,i,j,kfn,GRS)
             goto 9101 
          endif
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          if(knd==0)SK=gra(i,j)*1.d-5+1.5d0*gr*ksi(i,j)/r1
          if(knd==1)SK=gra(i,j)*1.d-5-0.5d0*gr*ksi(i,j)/r1
          rst=rst+SK/2.d0/pi*ds*(dwm(i,j)-dwm0)/dl**3
9101      continue
	  enddo
9100    continue
	enddo
	MolodeskyBLH=rst*1.d5
9002	return
      end
!--------------------------------------------------------------------------------
      real*8 function MolodeskySgn(BLH,gra,ksi,hgt,dwm,nlat,nlon,hd,knd,i0,j0,m,GRS)
      !细化核函数，计算BLH点的Hotine奇异积分
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
!-------------------------------------------------------------
      implicit none
	integer::knd,m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,gra(nlat,nlon),hgt(nlat,nlon),dwm(nlat,nlon),ksi(nlat,nlon)
	real*8::hd(6),pi,RAD,ds,mdr,r0,r1,r2,rst,rv,dwm0
	real*8::GRS(6),BLH(3),XYZ(3),rln(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,dl,SK,lon,lat,dg,dk,dw,gr,NFD(5)
!-----------------------------------------------------------------
      rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      dwm0=CGrdPntD2(BLH(2),BLH(1),dwm,nlat,nlon,hd)
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_RLAT(GRS,BLH,rln)
      rst=0.d0;mdr=rln(1)*rv*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),hgt,nlat,nlon,hd)
          dg=CGrdPntD2(BLH1(2),BLH1(1),gra,nlat,nlon,hd)
          dk=CGrdPntD2(BLH1(2),BLH1(1),ksi,nlat,nlon,hd)
          dw=CGrdPntD2(BLH1(2),BLH1(1),dwm,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
          dl=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          if(dl<mdr)goto 9101
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          if(knd==0)SK=dg*1.d-5+1.5d0*gr*dk/r1
          if(knd==1)SK=dg*1.d-5-0.5d0*gr*dk/r1
          rst=rst+SK/2.d0/pi*ds*(dw-dwm0)/dl**3
9101      continue
	  enddo
9100    continue
	enddo
	MolodeskySgn=rst
9002	return
      end
