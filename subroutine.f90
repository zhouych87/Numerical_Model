      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE bksub(jf,k1,k2,c)
      use commondat
      INTEGER:: jf,k1,k2
      REAL(16):: c(nci,ncj,nck) !Backsubstitution, used internally by solvde.
      INTEGER:: i,im,j,k,kp,nbf
      REAL(16):: xx
      nbf=ne-nb
      im=1
      do  k=k2,k1,-1   !Use recurrence relations to eliminate remaining dependences.
          if (k.eq.k1) im=nbf+1 !Special handling of rst point.
          kp=k+1
          do j=1,nbf
              xx=c(j,jf,kp)
              do i=im,ne
                  c(i,jf,k)=c(i,jf,k)-c(i,j,k)*xx
              end do 
          end do 
      end do 
      do k=k1,k2 !Reorder corrections to be in column 1.
          kp=k+1
          do  i=1,nb
              c(i,1,k)=c(i+nbf,jf,k)
          end do 
          do  i=1,nbf
              c(i+nb,1,k)=c(i,jf,kp)
          end do 
      end do 
      return
      END 
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE pinvs(ie1,ie2,je1,jsf,jc1,k,c,s)
      use commondat
      INTEGER:: ie1,ie2,jc1,je1,jsf,k
      REAL(16):: c(nci,ncj,nck),s(nsi,nsj)
      INTEGER,PARAMETER:: NMAX=10
      !Diagonalize the square subsection of the s matrix, and store the recursion coefficients in
      !c; used internally by solvde.
      INTEGER:: i,icoff,id,ipiv,irow,j,jcoff,je2,jp,jpiv,js1,indxr(NMAX)
      REAL(16):: big,dum,piv,pivinv,pscl(NMAX)
      je2=je1+ie2-ie1
      js1=je2+1
!      write(*,*) "Those row and column will be handled by pinvs"
!       write(*,*) "row S, and E,  column S,     B,      C-start point "
!       write(*,*) ie1,ie2,je1,je2
!       write(*,*) "at the front of pinve"
!       do  i=ie1,ie2; write(*,"(13(xg11.2))")  s(i,:) ;   end do
       
      do  i=ie1,ie2 !Implicit pivoting, as in x2.1.
          big=0.
          do  j=je1,je2
              if(abs(s(i,j)).gt.big) big=abs(s(i,j))
          end do 
          if(big.eq.0.) then 
             write(11,*) 'singular matrix, row all 0 in pinvs, row:',i 
             st=0 
             return
         end if 
          pscl(i)=1./big
          indxr(i)=0
!          write(*,*) i, indxr(i)
      end do 
!      write(*,*) "After compare"
!      do  i=1,6; write(*,"(13(xg11.2))")  s(i,:) ;   end do
      
      do  id=ie1,ie2
          piv=0.
          do  i=ie1,ie2 !Find pivot element.
              if(indxr(i).eq.0) then
                big=0.
                do  j=je1,je2
                  if(abs(s(i,j)).gt.big) then
                      jp=j
                      big=abs(s(i,j))
                  end if 
                end do 
                if(big*pscl(i).gt.piv) then
                  ipiv=i
                  jpiv=jp
                  piv=big*pscl(i)
                end if 
              end if 
          end do 
!           write(*,*) "After wwwwwww"
          if(s(ipiv,jpiv).eq.0.) then 
             write(11,*) 'singular matrix in pinvs'
             st=0 
             return
          end if 
          indxr(ipiv)=jpiv !In place reduction. Save column ordering.
          pivinv=1./s(ipiv,jpiv)
!          write(*,*) "find right ?",ipiv,jpiv,s(ipiv,jpiv)
          do  j=je1,jsf !Normalize pivot row.
              s(ipiv,j)=s(ipiv,j)*pivinv
          end do 
!          write(*,*) s(ipiv,jpiv)," =? 1"
          s(ipiv,jpiv)=1.
          do i=ie1,ie2 !Reduce nonpivot elements in column.
              if(indxr(i).ne.jpiv) then
                  if(s(i,jpiv).ne.0.) then
                      dum=s(i,jpiv)
                      do  j=je1,jsf
                          s(i,j)=s(i,j)-dum*s(ipiv,j)
                      end do 
                      s(i,jpiv)=0.
                  end if 
              end if 
          end do 
      end do 
            
      jcoff=jc1-js1 !js1=je2+1, jc1=1,set js1>1,and jsf>jsf-js1+1. Sort and store unreduced coefficients.
      icoff=ie1-je1 ! difference of row and column, is -3 
!      write(*,*) ie1,je1,jcoff , icoff
      do  i=ie1,ie2
          irow=indxr(i)+icoff ; 
          do  j=js1,jsf
!           write(*,*) "irow",i,irow,j,jcoff
              c(irow,j+jcoff,k)=s(i,j)
          end do 
      end do 
!      write(*,*) "After pivot"
!      do  i=1,ne; write(*,"(13(xg11.2))")  s(i,:) ;   end do
!      read(*,*) i 
      return
      END
      
      
      !!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,s)
      use commondat
      INTEGER:: ic1,iz1,iz2,jc1,jcf,jm1,jm2,jmf,jz1,jz2,kc 
      REAL(16):: c(nci,ncj,nck),s(nsi,nsj)
      !Reduce columns jz1-jz2 of the s matrix, using previous results as stored in the c matrix.
      !Only columns jm1-jm2,jmf are affected by the prior results. red is used internally by solvde.
      INTEGER:: i,ic,j,l,loff
      REAL(16):: vx
      loff=jc1-jm1
      ic=ic1
!      write(*,*) jz1,jz2, "columns in", iz1,iz2,"row need to be zeroed. "
!      write(*,*) "prior C,ic1,jz1,jz2",ic1,jm1,jm2,loff
      
!      write(*,*) "prior C matrix used by RED"
!      do  i=ic1,ic1+jz2-jz1; write(*,"(a,13(xg8.2))") "c",c(i,(jm1+loff):(jm2+loff),kc),c(i,jcf,kc); End do 
      
      do j=jz1,jz2 !Loop over columns to be zeroed.
          do  l=jm1,jm2 !Loop over columns altered.
              vx=c(ic,l+loff,kc)
              do  i=iz1,iz2 !Loop over rows.
                  s(i,l)=s(i,l)-s(i,j)*vx
              end do 
          end do 
          vx=c(ic,jcf,kc)
          do  i=iz1,iz2 !Plus final element.
              s(i,jmf)=s(i,jmf)-s(i,j)*vx
          end do 
          ic=ic+1
      end do 
!      do  i=1,ne; write(*,"(13(xg11.2))")  s(i,:) ;   end do
      return
      END
      
  