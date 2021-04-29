
      SUBROUTINE lnwab(ear, ne, param, ifl, photar, photer)

c      implicit none

      integer ne, ifl
      real ear(0:ne), param(2), photar(ne), photer(ne)
      real abs(ne), phot1(ne)

c    code to integrate over a log-normal distribution of column density
c
c    param(1) = nHmean (the mean column density) 10^22 cm^-2
c    param(2) = nHsigma (the dispersion of the column density) 10^22 cm^-2

c    local variables

      integer ie, ix, nx
      real nh, nhmean, nhsigma, x, dx, xmean
      real weight(31)

      nhmean=param(1)
      nhsigma=param(2)

c     If nhsigma = 0 then this is just the standard Wisconsin absorber

      IF ( nhsigma .EQ. 0 ) THEN
         CALL xsabsw(ear, ne, nhmean, ifl, photar, photer)
         RETURN
      ENDIF

c     set nh =1.e22 to calculate cold abosrption in this column density.
      nh=1.0 

c     check arrays are clear

      do ie=1,ne,1
         abs(ie)=0.0
         photar(ie)=0.0
         phot1(ie)=0.0
      end do

c     use the standard XSPEC wisconsin absorber but unlogged 
c     gives photar(E)=-nh*sigma(E)

      CALL xsabsw(ear, ne, nh, ifl, photar, photer)
c     use phot1(E)=-nh*sigma(E) instead of photar(E)=10**(-nh*sigma(E))

      phot1=log(photar)
c     set the total covering fraction unity

c     loop round values of Nh. lognh is initially set outside the loop
c     to avoid an undefined variable in the case when nnh=0.

c     use x = log(nh)/nhsigma as the variable
      dx = 0.2
      nx = 31 
c     1 + INT(6./dx)
      xmean = log(nhmean)/nhsigma
      x = xmean - 3.

c     calculate the Gaussian weight
      do ix = 1, nx
         weight(ix) = (erf((x+dx-xmean)/sqrt(2.))
     &      -erf((x-xmean)/sqrt(2.)))/2.
         x = x + dx
      end do

c     calculate the absorption
      x = xmean - 3.
      do ix = 1, nx
         nh=exp(nhsigma*(x+dx/2.0))
         do ie=1,ne,1
            abs(ie)=abs(ie)+exp(phot1(ie)*nh)*weight(ix)
         end do
         x = x + dx
      end do                    

      do ie=1,ne,1
         photar(ie)=abs(ie)
      end do

      return
      end
