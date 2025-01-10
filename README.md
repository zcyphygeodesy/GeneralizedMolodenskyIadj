## Fortran codes for Molodensky boundary value correction to Stokes on arbitrary shape surface
https://www.zcyphygeodesy.com/en/h-nd-151.html
## [Algorithm purpose]
&emsp;```From the gravity anomaly or gravity disturbance (mGal) model grid, height anomaly model grid, ellipsoidal height grid of the boundary surface and ellipsoidal height grid of the reference equipotential surface (for Stokes boundary problem), calculate the Molodensky I corrections of the gravity anomalies or gravity disturbances on the non-equipotential boundary surface on the ground or outside the Earth, thereby converting the Molodensky boundary value problem into a Stokes problem.```  
&emsp;```The boundary surface can be located at any altitude outside the geoid, and the shape of the boundary surface can be irregular.```
&emsp;```When the boundary surface is the ground surface and the reference equipotential surface is the geoid, the program calculates the classical Molodensky I (mGal).```
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg6_zltwYo2ZS0MjClDTjuCA.jpg)
## [Main program for test entrance]
    GeneralizedMolodenskyIadj.f90 
    Input parameters: knd - the type of the gravity field element to be corrected. Knd=0 for gravity anomaly, and knd=1 for gravity disturbance.
    Input parameters: dr - the integral radius (m).
    Input parameters: calcpntfl - the calculation point file name on the boundary surface. The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m), ......
    Input parameters: dbmhgrdfl - the ellipsoidal height grid file name of the boundary surface.
    Input parameters: dwmhgrdfl - the ellipsoidal height grid file name of the reference equipotential surface (for Stokes boundary problem).
    Input parameters: ksigrdfl - the height anomaly grid file name on the boundary surface.
    Input parameters: gravgrdfl - the gravity field element grid file name on the boundary surface. knd=0 for gravity anomaly, and knd=1 for gravity disturbance.
The same grid specifications required for the 4 input grid files above.
## (1) Algorithm module for numerical integral of Molodesky I item
    MolodeskyIntegral(calcpntfl,dbmhgrdfl,dwmhgrdfl,ksigrdfl,gravgrdfl,knd,dr)
    The output file reslt.txt, whose record format: Behind the record of the calculation point file, appends a column of the Molodensky boundary value correction.
## (2) Algorithm module for Molodesky boundary value correction to Stokes’
    real*8 function MolodeskyBLH(BLH,gra,ksi,hgt,dwm,nlat,nlon,hd,dr,knd,GRS)
    Input parameters BLH(3) - longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m) of the calculation point.
    Input parameters: hgt(nlat,nlon) - the ellipsoidal height grid (m) of the boundary surface.
    Input parameters: dwm(nlat,nlon) - the ellipsoidal height grid (m) of the reference equipotential surface.
    Input parameters: ksi(nlat,nlon) - the height anomaly grid (m) on the boundary surface.
    Input parameters: gra(nlat,nlon) - the gravity field element grid (mGal) of the boundary surface. knd=0 for gravity anomaly, and knd=1 for gravity disturbance.
    Input parameters: dr, hd(6) - the integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return: Molodesky boundary value correction (mGal).
## (3) Calculation module for the normal gravity field
    GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m²/s²), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (6) Other auxiliary modules
    PickRecord(line, kln, rec, nn); CGrdPntD2(lon,lat,dt,row,col,hd)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable file and all input and output data.
