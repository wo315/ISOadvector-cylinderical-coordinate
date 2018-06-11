    subroutine init
    use par 
    implicit none
    real(kind=dp)::Pooin,rho
    call grid
    !allocate variables and initicaltion
    allocate(U(0:NR+1,0:NZ+1),W(0:NR+1,0:NZ+1))
    allocate(Uo(0:NR+1,0:NZ+1),Wo(0:NR+1,0:NZ+1))
    allocate(Um(0:NR+1,0:NZ+1),Wm(0:NR+1,0:NZ+1))
    allocate(Nnormalcr(0:NR,0:NZ),Nnormalcz(0:NR,0:NZ))
    allocate(Po(1:NR,1:NZ),Pm(1:NR,1:NZ),P(1:NR,1:NZ))
    
    allocate(radv(1:NR,1:NZ),zadv(1:NR,1:NZ))
    allocate(radvo(1:NR,1:NZ),zadvo(1:NR,1:NZ))
    
    allocate(rdif(1:NR,1:NZ),zdif(1:NR,1:NZ))
    allocate(rdifo(1:NR,1:NZ),zdifo(1:NR,1:NZ))
    
    allocate(Rsource(1:NR,1:NZ),Zsource(1:NR,1:NZ),Uap(1:NR,1:NZ))
    allocate(Rsourceo(1:NR,1:NZ),Zsourceo(1:NR,1:NZ))
    
    allocate(Rcsf(1:NR,1:NZ),Zcsf(1:NR,1:NZ))

    allocate(Denc(1:NR,1:NZ),Visc(1:NR,1:NZ))
    !vof volume fraction
    allocate(Fnode(1:NR,1:NZ))
    allocate(F(-2:NR+2,-2:NZ+2))
    allocate(FN(-2:NR+2,-2:NZ+2),FO(-2:NR+2,-2:NZ+2))
    !vof 
    allocate(Cfl(1:NR,1:NZ),Nor(0:NR,0:NZ),Noz(0:NR,0:NZ))
    allocate(Coer(0:NR,0:NZ),Coez(0:NR,0:NZ))
    allocate(Coero(0:NR,0:NZ),Coezo(0:NR,0:NZ))
    !isoadvector
    allocate(fluxr(1:NR,1:NZ),fluxz(1:NR,1:NZ),vof_face_volr(1:NR,1:NZ),vof_face_volz(1:NR,1:NZ))
    allocate(check_bound(1:NR,1:NZ))
    !possion equation
    allocate(PAW(1:NR,1:NZ),PAE(1:NR,1:NZ))
    allocate(PAB(1:NR,1:NZ),PAT(1:NR,1:NZ),PAP(1:NR,1:NZ))
    allocate(PRHS(1:NR,1:NZ))

    !csf
    allocate(FAW(1:NR,1:NZ),FAE(1:NR,1:NZ))
    allocate(FAB(1:NR,1:NZ),FAT(1:NR,1:NZ),FAP(1:NR,1:NZ))
    allocate(Curvc(0:NR,0:NZ))
    !solver
    allocate(NAME(1:NR-1,1:NZ-1))

    !HF
    time=0._dp
    cyc=0
    dt=1e-5_dp
    nout=0
    !transform thetaA,thetaM
    thetaA=thetaA/180._dp*pi
    thetaM=thetaM/180._dp*pi
    
    !ÎïÐÔ²ÎÊýÓëËÙ¶È
    do i=0,NR+1
            do k=0,NZ+1
                Uo(i,k)=0.
                U(i,k)=0.
                Um(i,k)=0.
                Wo(i,k)=0.
                W(i,k)=0.
                Wm(i,k)=0.
            end do
    end do
     !ÎïÐÔ²ÎÊý,¸÷¸ö½Úµã
    do i=1,NR
            do k=1,NZ
                Denc(i,k)=0.
                Visc(i,k)=0.
            end do
    end do
    
    do i=1,NR
            do k=1,NZ
                Pm(i,k)=0.!ÐÞÕýÁ¿Îª0
                P(i,k)=0.
                Po(i,k)=0.
                radv(i,k)=0.
                zadv(i,k)=0.
                radvo(i,k)=0.
                zadvo(i,k)=0.
                Rdif(i,k)=0.
                Zdif(i,k)=0.
                Rdifo(i,k)=0.
                Zdifo(i,k)=0.
                Rsource(i,k)=0.
                Zsource(i,k)=0.
                Rsourceo(i,k)=0.
                Zsourceo(i,k)=0.
                Uap(i,k)=0.
                Rcsf(i,k)=0.
                Zcsf(i,k)=0.
          
                PAW(i,k)=0.
                PAE(i,k)=0.
                PAT(i,k)=0.
                PAB(i,k)=0.
                PAP(i,k)=0.
                PRHS(i,k)=0.
                
                !csf
                FAW(i,k)=0.
                FAE(i,k)=0.
                FAT(i,k)=0.
                FAB(i,k)=0.
                FAP(i,k)=0.
                Curvc(i,k)=0.
                !vof
                Cfl(i,k)=0.
                Fnode(i,k)=0.
            end do
    end do
    
     do i=1,NR
            do k=1,NZ
                PRHS(i,k)=0.
            end do
     end do
 
    !initial F
    do i=-2,NR+2
            do k=-2,NZ+2
                F(i,k)=0.
                FN(i,k)=F(i,k)
                FO(i,k)=F(i,k)
            end do
    end do

    do i=0,NR
            do k=0,NZ
                Nor(i,k)=0.
                Noz(i,k)=0.
                Coer(i,k)=0.
                Coez(i,k)=0.
                Coero(i,k)=0.
                Coezo(i,k)=0.
            end do
    end do


    if(set_compute==1) then
        call set_freesurface
    else
        call countinue
        call boundary_velocity
        call boundary_f
        call next 
    end if


    do k=1,NZ-1
            do i=1,NR-1
                number=(k-ZS)*(NR-1)+i-RS+1
                NAME(i,k)=number
            end do
    end do
    
  !get p*
    Pooin=0.
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    do i=RS,RE
           do k=ZS,ZE
            P(i,k)=P(i,k)-Pooin
           end do
    end do
    
    return
    end subroutine init

!method 1

subroutine set_freesurface
    use par 
    implicit none
    real::Heaviside
    real(kind=dp)::thick
    thick=sqrt(2.0)/2.0*min(dr,dz)
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                F(i,k)=(RC(i))**2+(ZC(k)-h0)**2
                !bubble
                F(i,k)=sqrt(F(i,k))-r0
                !water
!                 F(i,k)=r0-sqrt(F(i,k))

                F(i,k)=Heaviside(F(i,k),thick)
                  P(i,k)=epiron*2.0/r0*(1.0-F(i,k))
!                 P(i,k)=epiron*2.0/r0*F(i,k)
!                 if((RC(i))**2+(ZC(k)-h0)**2<r0**2) then
!                     W(i,k)=-1.58
!                     W(i,k+1)=-1.58
!                     P(i,k)=P(i,k)+Wdic*9.81*(sqrt(r0**2-(RC(i))**2)-ZC(k))
!                 end if
                !F(i,k)=max(min(1.0-(F(i,k)+thick)/2.0/thick,1.0),0.0)
            end do 
    end do 

    call boundary_f
    call boundary_velocity
    call next

    return

end subroutine  set_freesurface


!method two

subroutine set_freesurface2
    use par 
    implicit none
    real::Heaviside
    real(kind=dp)::thick
    integer::labelt(2),cellloc
    real(kind=dp)::Pointloc(4,2)
    thick=sqrt(2.0)/2.0*min(dr,dz)
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                labelt(1)=i
                labelt(2)=k
                call vertexloc(labelt,Pointloc)
                call isinterface(r0,h0,Pointloc,cellloc)
                if(cellloc==0) then
                    F(i,k)=(RC(i))**2+(ZC(k)-h0)**2
                    !bubble
                    !F(i,k)=sqrt(F(i,k))-r0
                    !water
                    F(i,k)=r0-sqrt(F(i,k))
                    
                    F(i,k)=Heaviside(F(i,k),thick)
                else if(cellloc==-1) then
                    F(i,k)=0.
                else if(cellloc==1) then
                    F(i,k)=1.
                end if 
                P(i,k)=epiron*2.0/r0*F(i,k)
                !F(i,k)=(RC(i))**2+(ZC(k)-h0)**2
                !F(i,k)=r0-sqrt(F(i,k))
                !F(i,k)=Heaviside(F(i,k),thick)
                !F(i,k)=max(min(1.0-(F(i,k)+thick)/2.0/thick,1.0),0.0)
            end do 
    end do 
    !do i=RS,RE
           ! do k=ZS,ZE
               ! if(F(i,k)>0.5) then
                   ! W(i,k)=-1.64
                   ! W(i,k+1)=-1.64
                !end if 
                !P(i,k)=epiron*2.0/r0*F(i,k)
            !end do 
    !end do 
    !do i=RS,RE
        !do k=ZS,ZE
           ! if(F(i,k)>0.5_dp) then
               ! P(i,k)=P(i,k)+Wdic*9.81*(sqrt(r0**2-RC(i)**2)+r0-Zc(k))
           ! end if 
        !end do 
   ! end do 

    call boundary_f
    call boundary_velocity
    call next

    return

end subroutine  set_freesurface2

real function   Heaviside(levvalue,thick)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp)::levvalue,thick
    real(kind=dp),parameter::pi=4_dp*atan(1.)
    if(levvalue<-thick) then
        Heaviside=0.
    else if(abs(levvalue)<=thick) then
        Heaviside=0.5*(1.0+levvalue/thick+1.0/pi*sin(pi*levvalue/thick))
    else
        Heaviside=1.0
    end if 
end function Heaviside


integer function newmod(inv,range)
    !help a array to behave like a linke list
    !test ok 2018_3-28
    implicit none
    integer::range,inv
    if(inv<=range) then
        newmod=inv
    else
        newmod=mod(inv,range)
    end if 
end function newmod 

subroutine isinterface(r1,h1,verloc,interfacestate)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::verloc(4,2)
    real(kind=dp),intent(in)::r1,h1
    integer::i,newmod
    logical::cutstate
    real(kind=dp)::distance(4),maxdistance,mindistance
    integer,intent(out)::interfacestate
    do i=1,4
        distance(i)=r1-sqrt(verloc(i,1)**2+(verloc(i,2)-h1)**2)
    end do 

    cutstate=.false.
    do i=1,4
        if(distance(i)*distance(newmod(i+1,4))<0.) then
            cutstate=.true.
        end if 
    end do 

    if(cutstate) then
        interfacestate=0
    else 
        maxdistance=max(distance(1),distance(2),distance(3),distance(4))
        mindistance=min(distance(1),distance(2),distance(3),distance(4))
        if(maxdistance<=0.) then
            interfacestate=-1
        end if 
        if(mindistance>=0.) then
            interfacestate=1
        end if 
    end if 

    return
end subroutine isinterface

 subroutine countinue
        use par
        implicit none
        open(11,file='save.data')
        do k=1,NZ-1
            do i=1,NR
                read(11,*)U(i,k)
            end do    
        end do

        do i=1,NR-1
            do k=1,NZ 
                read(11,*)W(i,k)
            end do 
        end do 

        do i=1,NR-1
            do k=1,NZ-1
                read(11,*)P(i,k),F(i,k)
            end do    
        end do
        close(11)
        return
    end subroutine countinue
