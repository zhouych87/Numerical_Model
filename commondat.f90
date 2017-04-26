module commondat
      implicit none
      INTEGER,PARAMETER:: M=301,ne=6,nb=3
      INTEGER,PARAMETER:: NCI=NE,NCJ=NE-NB+1,NCK=M+1
      INTEGER,PARAMETER:: NSI=NE,NSJ=2*NE+1,NYJ=NE,NYK=M
      REAL(16), public::  gtmp,h,slowc,conv,nc,nv,bgapv,V,DIELECT,kt,lti,ni
      REAL(16), DIMENSION(M):: X
      INTEGER, public:: mm,n,wout,werr,itmax,npoint,st,rmode
end module commondat

!!(mode,T,De,Dh,bgap,E0,Ed,t0,l)
