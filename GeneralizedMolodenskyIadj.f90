!  GeneralizedMolodenskyIadj.f90 
!
!  FUNCTIONS:
!  GeneralizedMolodenskyIadj - Entry point of console application.
!
!****************************************************************************
      program GeneralizedMolodenskyIadj
      implicit none
	character*800::calcpntfl,dbmhgrdfl,dwmhgrdfl,ksigrdfl,gravgrdfl
	integer knd
	real*8::dr
!---------------------------------------------------------------------
      knd=1!=0 for gravity anomaly, =1 for gravity disturbance
      !输入积分半径(m)
      dr=120.d3!Integral radius (m)
      !输入计算点文件名
      !Input the calculation point file on the boundary surface.
      write(calcpntfl,*)'dbmrga.txt'
      !输入边界面大地高格网文件名
      !Input the ellipsoidal height grid file of the boundary surface.
      write(dbmhgrdfl,*)'dbmhgt150s.dat'
      !输入参考等位面(Stokes边值归算面)大地高格网文件名
      !Input the ellipsoidal height grid file of the reference equipotential surface.
      write(dwmhgrdfl,*)'dwmhgt150s.dat'
      !输入边界面上高程异常格网文件名。
      !Input residual gravity field element grid file on the boundary  surface.
      write(ksigrdfl,*)'dbmGM1800150sksi.dat'
      !输入边界面上扰动场元格网文件名。mode=0 for gravity anomaly, mode=1 for gravity disturbance
      !Input anomalous gravity field element grid file on the boundary surface.
      write(gravgrdfl,*)'dbmGM1800150srga.dat'
      !程序要求4个输入格网文件有相同的格网规格
      !The same grid specifications required for the 4 input grid files above.
      write(*, *)"    Begin compulation......"
      call MolodeskyIntegral(calcpntfl,dbmhgrdfl,dwmhgrdfl,ksigrdfl,gravgrdfl,knd,dr)
      pause
      end
