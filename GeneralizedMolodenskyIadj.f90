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
      !������ְ뾶(m)
      dr=120.d3!Integral radius (m)
      !���������ļ���
      !Input the calculation point file on the boundary surface.
      write(calcpntfl,*)'dbmrga.txt'
      !����߽����ظ߸����ļ���
      !Input the ellipsoidal height grid file of the boundary surface.
      write(dbmhgrdfl,*)'dbmhgt150s.dat'
      !����ο���λ��(Stokes��ֵ������)��ظ߸����ļ���
      !Input the ellipsoidal height grid file of the reference equipotential surface.
      write(dwmhgrdfl,*)'dwmhgt150s.dat'
      !����߽����ϸ߳��쳣�����ļ�����
      !Input residual gravity field element grid file on the boundary  surface.
      write(ksigrdfl,*)'dbmGM1800150sksi.dat'
      !����߽������Ŷ���Ԫ�����ļ�����mode=0 for gravity anomaly, mode=1 for gravity disturbance
      !Input anomalous gravity field element grid file on the boundary surface.
      write(gravgrdfl,*)'dbmGM1800150srga.dat'
      !����Ҫ��4����������ļ�����ͬ�ĸ������
      !The same grid specifications required for the 4 input grid files above.
      write(*, *)"    Begin compulation......"
      call MolodeskyIntegral(calcpntfl,dbmhgrdfl,dwmhgrdfl,ksigrdfl,gravgrdfl,knd,dr)
      pause
      end
