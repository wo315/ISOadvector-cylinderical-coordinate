!boundary condition of contact_angle model:Numerical Analysis for Propellant Management in Rocket Tanks.pdf
!PhD-Numerical-prediction-of-two-fluid-systems-with-sharp-interfaces 
!Computational Fluid Dynamics of Dispersed Two-Phase Flows at
!High Phase Fractions------smooth_f_elliptic

!curvature compute:TWO-PHASE HEAT AND MASS TRANSFER MODELING: FLEXIBLE NUMERICAL METHODS FOR ENERGY
!ENGINEERING ANALYSES

subroutine csf
    use par 
    implicit none 

    real(kind=dp)::w1,W2,Cur,delta,H1,H2,CFc
    integer::count 

    CFc=0.98_dp

    do i=-1,NR+1  
        do k=-1,NZ+1
            FN(i,k)=F(i,k)
        end do
    end do
    
    call smooth_f

    call curvature_compute


    RS=2
    RE=NR-1
    ZS=1
    ZE=NZ-1
    
    Rcsfavge=0.
    count=0
    do i=RS,RE
            do k=ZS,ZE
                W1=sqrt(F(i,k)*(1._dp-F(i,k)))
                W2=sqrt(F(i-1,k)*(1._dp-F(i-1,k)))
                H1=1.0/(1.0-CFc)*(min(max(F(i,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                H2=1.0/(1.0-CFc)*(min(max(F(i-1,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                delta=(H1-H2)/dr

                Cur=(W1*Curvc(i,k)+W2*Curvc(i-1,k))/(W1+W2+1.e-20_dp)

                Rcsf(i,k)=epiron*delta*Cur
                if(abs(Rcsf(i,k))>1.e-12) then
                    Rcsfavge=Rcsfavge+Rcsf(i,k)
                    count=count+1
                end if 
            end do
    end do
    
    !get Rcsfavge
    Rcsfavge=Rcsfavge/float(count)

    RS=1
    RE=NR-1
    
    ZS=2
    ZE=NZ-1
    Zcsfavge=0.
    count=0
    do i=RS,RE
            do k=ZS,ZE
                W1=sqrt(F(i,k)*(1._dp-F(i,k)))
                W2=sqrt(F(i,k-1)*(1._dp-F(i,k-1)))
                H1=1.0/(1.0-CFc)*(min(max(F(i,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                H2=1.0/(1.0-CFc)*(min(max(F(i,k-1),CFc/2.),1.-CFc/2.)-CFc/2.)
                delta=(H1-H2)/dz

                Cur=(W1*Curvc(i,k)+W2*Curvc(i,k-1))/(W1+W2+1.e-20_dp)

                Zcsf(i,k)=epiron*delta*Cur
                if(abs(Zcsf(i,k))>1.0e-12) then
                    Zcsfavge=Zcsfavge+Zcsf(i,k)
                    count=count+1
                end if 
            end do
    end do

    Zcsfavge=Zcsfavge/float(count)
    
    return 
end subroutine csf 


subroutine curvature_compute
    use par
    implicit none
    real(kind=dp)::sgn
    real(kind=dp)::middle,top1,down1,top2,down2,top3,down3,Hzrr,Hzr
    integer::ii,kk,itin,itout
    real(kind=dp)::vp,Afw,Afe,Afb,Aft,Ffw,Ffe,Ffb,Fft,theta
    real(kind=dp)::Rcomponent,Zcomponent,mold,part11,part12,part21,part22
    real(kind=dp)::part2down,part21up,part22up,weight
    real(kind=dp)::sequence1(-1:1,-1:1),sequence2(-1:1,-1:1),sequence3(-1:1,-1:1)
    real(kind=dp),pointer,dimension(:,:)::Nnormalfr,Nnormalfz
    real(kind=dp),pointer,dimension(:,:,:,:)::smoothvalue2,smoothvalue1
    allocate(Nnormalfr(0:NR,0:NZ),Nnormalfz(0:NR,0:NZ))
    allocate(smoothvalue1(0:NR,0:NZ,0:2,0:2),smoothvalue2(0:NR,0:NZ,0:2,0:2))
    !nodes F value
    !initial

    do i=0,NR
        do k=0,NZ
            Nnormalcr(i,k)=0.
            Nnormalcz(i,k)=0.
            Nnormalfr(i,k)=0.
            Nnormalfz(i,k)=0.
        end do 
    end do 

    do ii=0,2
        do kk=0,2
            do i=0,NR
                do k=0,NZ
                    smoothvalue1(i,k,ii,kk)=0.
                    smoothvalue2(i,k,ii,kk)=0.
                end do 
            end do 
        end do      
    end do
     !...................................................................
     !step1...........compute the center unit normal vector with CNC
 	do i=1,NR 
 		do k=1,NZ 
 			Fnode(i,k)=0.25*(FN(i,k)+FN(i-1,k)+FN(i,k-1)+FN(i-1,k-1))
 		end do 
 	end do 

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

 	do i=RS,RE 
 		do k=ZS,ZE

                vp=dr*dz
                Afw=dz
                Afe=dz
                Afb=dr
                Aft=dr
     

 			Ffw=0.5_dp*(Fnode(i,k)+Fnode(i,k+1))
 			Ffe=0.5_dp*(Fnode(i+1,k)+Fnode(i+1,k+1))
 			Ffb=0.5_dp*(Fnode(i,k)+Fnode(i+1,k))
 			Fft=0.5_dp*(Fnode(i,k+1)+Fnode(i+1,k+1))
            Rcomponent=(Ffe-Ffw)/dr
            Zcomponent=(FFt-Ffb)/dz

 			mold=sqrt(Rcomponent**2+Zcomponent**2)
 			Nnormalcr(i,k)=Rcomponent/(mold+1.0e-20)
 			Nnormalcz(i,k)=Zcomponent/(mold+1.0e-20)

 		end do 
 	end do 
    !...................................................................
    !setp2.......smooth the center unit normal
    !extend the center value
    if(BCXLEF==3) then
        call sy_extend(Nnormalcr)
    else
        call center_extend(Nnormalcr)
    end if 
    call center_extend(Nnormalcz)

    !no we get the smooth unit normal vector
    !...................................................................
    !step1 get the face unit normal component
    !use normal center vector to compute the normal face vector
    do i=1,NR 
        do k=1,NZ-1 
           Rcomponent=0.5_dp*(Nnormalcr(i,k)+Nnormalcr(i-1,k))
           Zcomponent=0.5_dp*(Nnormalcz(i,k)+Nnormalcz(i-1,k))
           mold=sqrt(Rcomponent**2+Zcomponent**2)
           Nnormalfr(i,k)=Rcomponent/(mold+1.0e-20) 
           !store the face's normalr
        end do 
    end do 

    do k=1,NZ  
        do i=1,NR-1
           Rcomponent=0.5_dp*(Nnormalcr(i,k)+Nnormalcr(i,k-1))
           Zcomponent=0.5_dp*(Nnormalcz(i,k)+Nnormalcz(i,k-1))
           mold=sqrt(Rcomponent**2+Zcomponent**2)
           Nnormalfz(i,k)=Zcomponent/(mold+1.0e-20)
           !store the face's normalz
        end do 
    end do 
    !...................................................................
    !step 2 compute the center curvature
    !use the face normal value to compute the center curvature

    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1  

    do i=RS,RE
        do k=ZS,ZE
                vp=dr*dz
                Afw=dz
                Afe=dz
                Afb=dr
                Aft=dr

            Ffw=Nnormalfr(i,k)
            Ffe=Nnormalfr(i+1,k)
                
            Ffb=Nnormalfz(i,k)
            Fft=Nnormalfz(i,k+1)
            !boundary condition
            Curvc(i,k)=-(Afe*ffe-Afw*ffw+Aft*Fft-Afb*Ffb)/vp
        end do 
    end do 

    !...................................................................
    !smooth the center curvature
    !step1 extend the center curvature 
    call center_extend(Curvc) 
    !step 2 smooth the center curvature use the CNC_scheme scheme
    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do itout=0,1
        do itin=0,1

            do i=0,NR
                do k=0,NZ
                    smoothvalue1(i,k,0,0)=Curvc(i,k) !initial value
                end do
            end do 

            do i=RS,RE
                do k=ZS,ZE
                    if(itout==0) then
                        part11=smoothvalue1(i,k,0,0)
                    else
                        part11=smoothvalue1(i,k,itin,itout-1)
                    end if 
                    weight=(4.0*F(i,k)*(1.0-F(i,k))+1.0e-12)**0.25
                    do ii=-1,1
                        do kk=-1,1
                            sequence1(ii,kk)=(4.0*F(i+ii,k+kk)*(1.0-F(i+ii,k+kk))+1.0e-12)**0.25
                            sequence2(ii,kk)=sequence1(ii,kk)*smoothvalue1(i+ii,k+kk,itin,itout)
                        end do 
                    end do 
                    call CNC_scheme(sequence1,part2down)
                    call CNC_scheme(sequence2,part21up)
                    part21=part21up/part2down
                    smoothvalue1(i,k,itin+1,itout)=weight*part11+(1.0-weight)*part21
                end do 
            end do 
            !extend to get smooth values
            do k=ZS,ZE
                i=RS-1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i+1,k,itin+1,itout)
                i=RE+1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i-1,k,itin+1,itout)
            end do 

            do i=RS-1,RE+1
                k=ZS-1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i,k+1,itin+1,itout)
                k=ZE+1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i,k-1,itin+1,itout)
            end do 
            !update values
            do i=0,NR
                do k=0,NZ
                   Curvc(i,k)=smoothvalue1(i,k,itin+1,itout)
                end do
            end do 
            !update curvc complete
      

        end do ! end of inner loop
    end do ! end of outer loop

    
    deallocate(Nnormalfr,Nnormalfz)!Nnormalcr,Nnormalcz
    deallocate(smoothvalue2,smoothvalue1)
    return

end subroutine curvature_compute

        
subroutine center_extend(center)
	use par
    implicit none
    real(kind=dp)::center(0:NR,0:NZ)

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do k=ZS,ZE
    	i=RS-1
    	center(i,k)=center(i+1,k)
    	i=RE+1
    	center(i,k)=center(i-1,k)
    end do 

    do i=RS-1,RE+1
    	k=ZS-1
    	center(i,k)=center(i,k+1)
    	k=ZE+1
    	center(i,k)=center(i,k-1)
    end do 
    return

end subroutine center_extend 

subroutine sy_extend(center)
    use par
    implicit none
    real(kind=dp)::center(0:NR,0:NZ)

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do k=ZS,ZE
        i=RS-1
        center(i,k)=-center(i+1,k)
        i=RE+1
        center(i,k)=center(i-1,k)
    end do 

    do i=RS-1,RE+1
        k=ZS-1
        center(i,k)=center(i,k+1)
        k=ZE+1
        center(i,k)=center(i,k-1)
    end do 
    return

end subroutine sy_extend 

subroutine CNC_scheme(sequence,centervalue)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::sequence(-1:1,-1:1)
    real(kind=dp),intent(out)::centervalue

    centervalue=4.*sequence(0,0)+sequence(-1,1)+sequence(-1,-1)+sequence(1,-1)+&
    sequence(1,1)+2.0*(sequence(-1,0)+sequence(1,0)+sequence(0,-1)+sequence(0,1))

    centervalue=centervalue/16.

    return 
end subroutine CNC_scheme


subroutine smooth_f_csf
    use par 
    implicit none 
    integer::it,item,ii,kk
    real(kind=dp),pointer,dimension(:,:,:)::smoothvalue
    real(kind=dp)::part1,part2,CKS,sequence(-1:1,-1:1)
    allocate(smoothvalue(0:NR,0:NZ,0:2))

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do i=0,NR
        do k=0,NZ
            do ii=0,2
                smoothvalue(i,k,ii)=0.
            end do 
        end do 
    end do 

    it=0
    item=2

    do while(it<item)
        it=it+1
        !initial value
        do i=0,NR
            do k=0,NZ
                    smoothvalue(i,k,0)=FN(i,k) 
            end do
        end do 

        do i=RS,RE
            do k=ZS,ZE

                part1=smoothvalue(i,k,it)

                CKS=0.5

                do ii=-1,1
                    do kk=-1,1
                        sequence(ii,kk)=smoothvalue(i+ii,k+kk,it)
                    end do 
                end do 

                call CNC_scheme(sequence,part2)

                smoothvalue(i,k,it+1)=CKS*part1+(1.0-CKS)*part2

                end do 
        end do 
            !extend to get smooth values
        do k=ZS,ZE
                i=RS-1
                smoothvalue(i,k,it+1)=smoothvalue(i+1,k,it+1)
                i=RE+1
                smoothvalue(i,k,it+1)=smoothvalue(i-1,k,it+1)
        end do 

        do i=RS-1,RE+1
                k=ZS-1
                smoothvalue(i,k,it+1)=smoothvalue(i,k+1,it+1)
                k=ZE+1
                smoothvalue(i,k,it+1)=smoothvalue(i,k-1,it+1)
        end do 
            !update values
        do i=0,NR
            do k=0,NZ
                   FN(i,k)=smoothvalue(i,k,it+1)
            end do
        end do 
            !update curvc complete
      

    end do ! end of so while 
    deallocate(smoothvalue)
    return
end subroutine smooth_f_csf






