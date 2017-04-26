      program main
      use commondat
      ! USES plgndr,solvde
      !Sample program using solvde. Computes eigenvalues of spheroidal harmonics Smn(x; c)
      !for m   0 and n   m. In the program, m is mm, c2 is c2, 
      implicit none
      include 'mpif.h'
      INTEGER:: is,k,indexv(NE),mode,emode
      REAL(16):: deriv,q1,c(NCI,NCJ,NCK),s(NSI,NSJ),ac 
      REAL(16):: scalv(NE),y(NE,M),g,T,De,Dh,bgap,E0,Ed,LT,D
      integer :: myid,ierr,npcs,status(MPI_STATUS_SIZE) 
      real :: start, finish
!       T=298 kt=25.6meV
!      T0=57  = 273.15+57=330.15K
      call cpu_time(start)
      call init() 
      open(11,file="error.out",status="replace") 
      !!
        open(301,file="in",status="old") 
        read(301,*) mode ,emode, DIELECT ,rmode! 1: LT;  2:LT-g; 3 rnp, 4: rnp-g,
        write(*,*) mode 
        if (mode <=3) then 
           write(*,*) "Light intensity is uniform in the cell in this simulation"
        else if (mode>=4 .and. mod(mode,2)==0) then 
           write(*,*) "Light intensity is exponential and the light come in from TiO2 side"
        else if (mode>=4 .and. mod(mode,2)==1) then 
           write(*,*) "Light intensity is exponential and the light come in from Spiro-OMeTAD side"
        end if 
        if (mode ==1 .or. mode==4 .or. mode==5) then 
           write(*,*) "Recombination use light time (SRH) model"
        else if (mode ==2 .or. mode==6 .or. mode==7) then 
           write(*,*) "Direct recombination, r calculated from uniformly distributed carriers."
        else if (mode ==3 .or. mode==8 .or. mode==9) then 
           write(*,*) "Direct recombination, r calculated from exponentially distributed carriers."
        end if 
       ! mode>10: diffusion length change with temperature change 
        read(301,*) itmax,npoint ! maximun iteration, discreting points 
        read(301,*) D     ! thickness of solar cell
             D=D*1.0E-7 ! from nm to cm 
        read(301,*) de,dh ! diffusion coefficient of electron and hole
        read(301,*) bgap,LT ! bgap for absoprtion (/cm^2) and lifetime(ns)
        read(301,*) ac  ! absortption length cm 
        read(301,*) nc,nv  ! a
        read(301,*) bgapv ! ! bgap for output voltage 
        read(301,*) slowc,conv 
        read(301,*) e0,ed  ! field on two sides
        read(301,*)  T   ! working temperature 
        read(301,*) lti ! light intensity 
!        read(301,*)  
        close(301) 
             kt=0.0256*t/298 
        if (mode>10) then 
           de=de*30.15/abs(330.15-t)
           dh=dh*30.15/abs(330.15-t) 
           mode=mode-10 
        end if 
 !       write(*,*) de,dh
        call IV(mode,T,De,Dh,bgap,E0,Ed,LT,D,ac,Emode)
!        write(*,*) e0,ed
!        read(*,*)
      call MPI_Finalize(ierr) 
      call cpu_time(finish) 
      print '("Time = ",f15.3," seconds.")',finish-start
     close(11)
      end program main 


      SUBROUTINE IV(mode,T,De,Dh,bgap,E0,Ed,LT,D,ac,Emode)
      use commondat
      ! USES plgndr,solvde
      !Sample program using solvde. Computes eigenvalues of spheroidal harmonics Smn(x; c)
      !for m   0 and n   m. In the program, m is mm, c2 is c2, 
      implicit none
      include 'mpif.h'
      INTEGER:: is,k,indexv(NE),mode ,emode
      REAL(16):: c(NCI,NCJ,NCK),s(NSI,NSJ),ac
      REAL(16):: D,scalv(NE),y(NE,M),T,De,Dh,bgap,E0,Ed,LT,tmp
      integer :: myid,ierr,npcs,status(MPI_STATUS_SIZE) 

!      n0=3.97E+18*1.6E-19 !Nc  integral of DOS
!       n0=3.97E+18*1.6E-19
!       ni=n0*exp(-bgap/(2*kt))
       ni=1.0E+14*1.6E-19
!      write(*,*) "please notice that the DOS of perovskite is 10^21 "
!      write(*,"(a,2xg15.6)") "ni;", ni
      gtmp=-2.2E+17*bgap+5.115E+17
      gtmp=lti*gtmp
!       write(*,"(2(xg15.7))") (exp(1.0)-1)**2/(2*LT*LT*1.6E-19*gtmp/D), 3.18761/(1.6E-19*5.7E+4*gtmp*Lt)
            
      if (mode == 2 .or. mode==6 .or. mode==7) then 
          write(*,"(a,g15.7)") " rnp recombination coefficient", Lt !LT=1.03E-9/1.6E-19 !4.62E-9 FOR 380  LT=(exp(1.0)-1)**2/(LT*LT*1.6E-19*gtmp/D) !!rnp 
          LT=6.67E-04/LT**2
          LT=lt/1.6E-19
      else if (mode == 3 .or. mode==8 .or. mode==9) then 
         write(*,"(a,g15.7)") " rnp recombination coefficient", Lt !LT=4.62E-9/1.6E-19 !lt=3.18761/(1.6E-19*5.7E+4*gtmp*Lt)             !!!!!rnp2
          LT=6.67E-04/LT**2
          LT=lt/1.6E-19
      else 
            LT=LT*1.0E-9    ! t0
         write(*,"(a,g15.7)") "t0 model lifetime", Lt
    end if 

      h=D/(M-1) ! discrete 
      indexv(1)=1
      indexv(2)=2
      indexv(3)=3
      indexv(4)=4
      indexv(5)=5
      indexv(6)=6

     tmp=bgap/D

      call MPI_INIT( ierr )     
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )     
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npcs, ierr )
       if (myid== 0)  write (*,"(a)") "   it    V(mV)    J(mA)        n         p       E" 
       
       do is=myid,npoint,npcs 
!       do is=npoint-myid,1,-npcs 
            st =1
            y=0;s=0;c=0
!            write(*,*) is,'begin st',st
             v=1.3*is/npoint - bgap ! v=0-1300,v=
             if (emode==1) then 
               E0=-0.05*(V+bgap)/D!+650!  -0.05*1.0/D;  !! V/cm
               Ed=-0.05*(V+bgap)/D!+650!+0.05*0.3*tmp! -0.05*1.0/D;  !! V/cm
!               write(*,*) "Please notice the code of this part!!"
!               Ed=V/D
!               E0=ED
!              write(*,"(2(xg15.5))",advance="no") (v+bgap)*1000, e0 
            end if 
            do k=1,M !Initial guess.
                 x(k)=(k-1)*h
                 y(6,k)=(D-x(k))*0.01/D  ! jn
                 y(2,k)=0.01*x(k)/D       !jp
                 y(3,k)=E0
                 y(4,k)=n0*exp(v/(2*kt))*EXP(x(k)/D) !n 
                 y(5,k)=n0*exp(v/(2*kt))*EXP(1-x(k)/D) !p
!                 y(6,k)=0.3
!                 y(7,k)=0.01
                 y(1,k)=n0*exp(v/(2*kt))
             end do 
!             write(*,*) v
             scalv(6)=0.01
             scalv(2)=0.01
             scalv(3)=E0
             scalv(4)=n0*exp(v/(2*kt))
             scalv(5)=n0*exp(v/(2*kt))
!             scalv(6)=1
!             scalv(7)=0.01
             scalv(1)=n0*exp(v/(2*kt)) 
!             read(*,*) k 
             call solvde(scalv,indexv,y,c,s,mode,T,De,Dh,bgap,E0,Ed,LT,D,ac)
!             write(*,*) is,'end st:' ,st
      end do 
      return
      END
    
          
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE deq(k,k1,k2,jsf,is1,isf,indexv,s,y,mode,T,De,Dh,bgap,E0,Ed,LT,D,ac)
 use commondat
 implicit none
 INTEGER:: is1,isf,jsf,k,k1,k2,indexv(nyj),i,mode
 REAL(16):: s(nsi,nsj),y(nyj,nyk),tmp1,tmp2,tmp,ac,D,g,T,De,Dh,bgap,E0,Ed,LT
!      REAL(16):: temp,temp2 
!      ac=5.7E+4
      D=h*(M-1)
      s=0 
!      write(*,*) "k,  jsf , y31,y41"
!      write(*,*) k,  jsf,y(1,1),y(2,1)
   
  if(k.eq.k1) then !Boundary condition at x=0, TiO2
              s(4,ne+3)=1.0
              s(5,ne+2)=1.0;      
              s(6,ne+4)=1.0;        s(6,ne+1)=-1.0  ! n0 
              
              s(4,jsf)=y(3,1)-E0 
              s(5,jsf)=y(2,1) !jp =0 
              s(6,jsf)=y(4,1)-y(1,1) !n-n0 
  else if(k.gt.k2) then !Boundary conditions at x=M, Spiro.
              s(1,ne+3)=1.0 !Equation (17.4.32).
              s(2,ne+6)=1.0;  
              s(3,ne+5)=1.0;   s(3,ne+1)=exp(v/(kt))*n0*n0/(y(1,m)*y(1,m))
              
              s(1,jsf)=y(3,m)-Ed
              s(2,jsf)=y(6,m) !jn=0
              s(3,jsf)=y(5,m)-exp(v/(kt))*n0*n0/y(1,m)
!              do  i=1,6
!               write(*,"(13(g11.4))")  s(i,:) 
!              end do
  else                 !Interior point.
     if (mode<=3 ) then                                      ! uniform generation
          g=1.6E-19*gtmp/D 
      else if (mode==5 .or. mode==7 .or. mode==9 ) then  !exponential generation Light form spiro x=d
          g=1.6E-19*gtmp*ac*exp(-ac*(M-k+0.5)*h)
      else if (mode==4 .or. mode==6 .or. mode==8 ) then  !exponential generation Light form TiO2 x=0
          g=1.6E-19*gtmp*ac*exp(-ac*(k-0.5)*h)
      else 
         write(*,*) "mode error mode should between 1-9.  Current mode is", mode 
         stop 
      end if 
      
        tmp1=y(4,k)+y(4,k-1);tmp2=y(5,k)+y(5,k-1);tmp=tmp1+tmp2
      if (mode==1 .or. mode==4 .or. mode==5) then                        !! t0 model
          s(6,4)=-0.5*h*(tmp2/tmp-(tmp1*tmp2)/tmp**2)/LT 
          s(6,ne+4)= s(1,4)
          s(6,5)=-0.5*h*(tmp1/tmp-(tmp1*tmp2)/tmp**2)/LT 
          s(6,ne+5)=s(1,5)
          s(6,6)=-1.0
          s(6,ne+6)=1.0
          s(6,jsf)=y(6,k)-y(6,k-1)-h*(tmp1*tmp2/((tmp1+tmp2)*2*LT)-g)
          
          s(2,4)=-s(6,4)
          s(2,ne+4)=s(2,4)
          s(2,5)=-s(6,5)
          s(2,ne+5)=s(2,5)
          s(2,2)=-1.0
          s(2,ne+2)=1.0
          s(2,jsf)=y(2,k)-y(2,k-1)+h*(tmp1*tmp2/((tmp1+tmp2)*2*LT)-g)
      else                             !!! rnp and rnp2 model , the difference is the r coefficient
          s(6,6)=-1.0
          s(6,ne+6)=1.0
          s(6,4)=-0.25*h*tmp2*LT 
          s(6,ne+4)= s(6,4)
          s(6,5)=-0.25*h*tmp1*LT
          s(6,ne+5)=s(6,5)
          if(rmode==1) s(6,jsf)=y(6,k)-y(6,k-1)-h*((0.25*tmp1*tmp2*LT-ni*ni)-g)  ! r(n*p-ni*ni)
          if(rmode==0) s(6,jsf)=y(6,k)-y(6,k-1)-h*((0.25*tmp1*tmp2*LT)-g)  ! rnp
          
          s(2,4)=0.25*h*tmp2*LT
          s(2,ne+4)=s(2,4)
          s(2,5)=0.25*h*tmp1*LT
          s(2,ne+5)=s(2,5)
          s(2,2)=-1.0
          s(2,ne+2)=1.0
          if(rmode==1) s(2,jsf)=y(2,k)-y(2,k-1)+h*((0.25*tmp1*tmp2*LT-ni*ni)-g)
          if(rmode==0 ) s(2,jsf)=y(2,k)-y(2,k-1)+h*((0.25*tmp1*tmp2*LT)-g)
       end if 
          s(3,4)=0.5*h/(DIELECT*8.854E-14) ! tmp1=y4+y4 cm 
          s(3,ne+4)=0.5*h/(DIELECT*8.854E-14)
          s(3,5)=-0.5*h/(DIELECT*8.854E-14)
          s(3,ne+5)=-0.5*h/(DIELECT*8.854E-14)
          s(3,3)=-1.0
          s(3,ne+3)=1.0 
          s(3,jsf)=y(3,k)-y(3,k-1)-0.5*h*(tmp2-tmp1)/(DIELECT*8.854E-14)
!         s(3,jsf)=y(3,k)-y(3,k-1)+0.5*h*(tmp2-tmp1)/(DIELECT*8.854E-14)
          
          s(4,4)=-1.+h*(y(3,k)+y(3,k-1))/(4*kt) 
          s(4,ne+4)=1.+h*(y(3,k)+y(3,k-1))/(4*kt) 
          s(4,6)=-.5*h/de
          s(4,ne+6)=-.5*h/de
          s(4,3)=h*tmp1/(4*kt)
          s(4,ne+3)=h*tmp1/(4*kt)
          s(4,jsf)=y(4,k)-y(4,k-1)-0.5*h*((y(6,k)+y(6,k-1))/de-tmp1*(y(3,k)+y(3,k-1))/(2*kt))
          
          s(5,5)=-1.-h*(y(3,k)+y(3,k-1))/(4*kt) 
          s(5,ne+5)=1.-h*(y(3,k)+y(3,k-1))/(4*kt) 
          s(5,2)=0.5*h/dh
          s(5,ne+2)=0.5*h/dh
          s(5,3)=-h*tmp2/(4*kt)
          s(5,ne+3)=-h*tmp2/(4*kt) 
          s(5,jsf)=y(5,k)-y(5,k-1)+0.5*h*((y(2,k)+y(2,k-1))/dh-tmp2*(y(3,k)+y(3,k-1))/(2*kt)) 
          
          s(1,1)=-1; s(1,ne+1)=1  
          s(1,jsf)=y(1,k)-y(1,k-1)
     end if 
!          do  i=1,6
!               write(*,"(13(g11.4))")  s(i,:) 
!         end do
      return
      END
      
     SUBROUTINE init() 
     WRITE(*,*) "This code is written by Yecheng Zhou."
     WRITE(*,*) "Please check information and updates on https://github.com/zhouych87/Numerical_Model"
     WRITE(*,*) "Details are shown in Phys. Chem. Chem. Phys., 18(6), 4476–4486. "
     WRITE(*,*) "Please cite it when you using my code"
     WRITE(*,*) "Any question please send email to zhouych87@gmail.com"
     END 
      
      SUBROUTINE solvde(scalv,indexv,y,c,s,mode,T,De,Dh,bgap,E0,Ed,LT,D)
       use commondat
      implicit none
      INTEGER:: indexv(nyj),mode
      REAL(16):: c(nci,ncj,nck),s(nsi,nsj),scalv(nyj),y(nyj,nyk),T,De,Dh,bgap,E0,Ed,LT,D,cp
      INTEGER:: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9,jc1,jcf,jv,k,k1,k2,km,kp,nvars,kmax(ne)
      REAL(16):: err,errj,fac,vmax,vz,ermax(ne)
      character(20)::flname 
      k1=1 ! Set up row and column markers.
      k2=m
      nvars=ne*m
      j1=1
      j2=nb
      j3=nb+1
      j4=ne
      j5=j4+j1
      j6=j4+j2
      j7=j4+j3
      j8=j4+j4
      j9=j8+j1
      ic1=1
      ic2=ne-nb
      ic3=ic2+1
      ic4=ne
      jc1=1
      jcf=ic3
      do it=1,itmax !Primary iteration loop.
          k=k1 !Boundary conditions at first point.
          call deq(k,k1,k2,j9,ic3,ic4,indexv,s,y,mode,T,De,Dh,bgap,E0,Ed,LT,D)
          call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s) 
          if (st==0) then ; write(11,'(a,xf6.2)') "head voltage", bgap+V; return; end if  
          do  k=k1+1,k2 ! Finite difference equations at all point pairs.
            kp=k-1
!            write(*,*) "points:", k  
            call deq(k,k1,k2,j9,ic1,ic4,indexv,s,y,mode,T,De,Dh,bgap,E0,Ed,LT,D)
            call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp, c,s)
            call pinvs(ic1,ic4,j3,j9,jc1,k,c,s)
          if (st==0) then ; write(11,'(a,xf6.2)') "mid voltage", bgap+V; return; end if  
          end do 
          k=k2+1 !Final boundary conditions.
          call deq(k,k1,k2,j9,ic1,ic2,indexv,s,y,mode,T,De,Dh,bgap,E0,Ed,LT,D)
          call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s)
          call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
          if (st==0) then ; write(11,'(a,xf6.2)') "tail voltage", bgap+V; return; end if  
          call bksub(jcf,k1,k2,c) !Backsubstitution.
          err=0.
          do  j=1,ne !Convergence check, accumulate average error.
              jv=indexv(j)
              errj=0.
              km=0
              vmax=0.
              do  k=k1,k2 !Find point with largest error, for each dependent variable.
                  vz=abs(c(jv,1,k))
                  if(vz.gt.vmax) then
                  vmax=vz
                  km=k
                  end if 
                  errj=errj+vz
              end do 
              scalv(j)=maxval(abs(y(j,:)))
if (npoint==1)              write(*,"(a,xi0,2(xg15.4),xi0,xg15.4)") "error for equation:",j,errj,errj/scalv(j),km,vmax/scalv(j)
              err=err+errj/scalv(j) !Note weighting for each dependent variable.
              ermax(j)=c(jv,1,km)/scalv(j)
              kmax(j)=km 
!if (npoint==1)         write(*,'(i0,xg13.6)')    km,ermax(j)
          end do  
           err=err/nvars
          fac=slowc/max(slowc,err) !Reduce correction applied when error is large.
          if (npoint==1)          write(*,"(i6,2(xg13.6))") it,err,fac!Summary of corrections for this step. Point with largest error 
      !    for each variable can be monitored by writing out kmax and ermax.

          if(err.lt.conv .and. y(4,1)> 0 .and. y(5,m)>0) then 
!             write(flname,"(a,i0,a)") 'v', NINT((bgap*1000+25.6*t/298*log(y(4,1)*y(5,m)/n0**2))*10),'.out'
!             open(111,file=trim(flname),status="replace") 
!             write(111,"(a)") "  x     connect-n      Jp       E        N      P     Jn        Vacc       Ec       El"
!             cp=0 
!             do  k=1,m 
!                  cp=cp+h*y(3,k)
!         write(111,"(10(xg15.8))") x(k),(y(j,k),j=1,ne),cp,-3.94+kt*log(y(4,k)/n0),-5.50-kt*log(y(5,k)/n0)
!             end do 
!             close(111)
           write (*,"(i6, 2(xf15.8),3(xg15.8))") it, bgap*1000+25.6*t/298*log(y(4,1)*y(5,m)/n0**2),y(6,1)*1000,y(4,1),y(5,m),y(1,1)
             return
          end if 
          
          do  j=1,ne !Apply corrections.
              jv=indexv(j)
              do  k=k1,k2
                  y(j,k)=y(j,k)-fac*c(jv,1,k)
!                 if (j==4 .or. j==5 ) y(j,k)=abs(y(j,k))
              end do 
          end do 

      end do 
      write(11,*) 'itmax exceeded in solvde' !Convergence failed.
!      stop 
      return
      END 
      
