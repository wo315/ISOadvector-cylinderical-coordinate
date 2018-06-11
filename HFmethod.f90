subroutine csf
    use par 
    implicit none 

    real(kind=dp)::w1,W2,Cur,delta,H1,H2,CFc
    integer::count 

    CFc=0.98_dp

    !do i=-1,NR+1  
        !do k=-1,NZ+1
            !FN(i,k)=F(i,k)
        !end do
    !end do
    
    !call smooth_f

    call curvature_compute


    RS=2
    RE=NR-1
    ZS=1
    ZE=NZ-1
    
    Rcsfavge=0.
    count=0
    do i=RS,RE
            do k=ZS,ZE
                W1=(4.0*F(i,k)*(1._dp-F(i,k)))**0.25
                W2=(4.0*F(i-1,k)*(1._dp-F(i-1,k)))**0.25
                ! H1=1.0/(1.0-CFc)*(min(max(F(i,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                !H2=1.0/(1.0-CFc)*(min(max(F(i-1,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                if(F(i,k)>0.5) then
                    H1=1.0
                else
                    H1=0.0
                end if 
                if(F(i-1,k)>0.5) then
                    H2=1.0
                else
                    H2=0.
                end if 
                delta=(H1-H2)/dr

                Cur=(W1*Curvc(i,k)+W2*Curvc(i-1,k))/(W1+W2+1.e-20_dp)

                if((2.0*F(i,k)-1.0)*(2.0*F(i-1,k)-1.0)<0.) then
                    Rcsf(i,k)=epiron*delta*Cur
                else
                    Rcsf(i,k)=0.
                end if 

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
                W1=(4.0*F(i,k)*(1._dp-F(i,k)))**0.25
                W2=(4.0*F(i,k-1)*(1._dp-F(i,k-1)))**0.25
               ! H1=1.0/(1.0-CFc)*(min(max(F(i,k),CFc/2.),1.-CFc/2.)-CFc/2.)
                !H2=1.0/(1.0-CFc)*(min(max(F(i,k-1),CFc/2.),1.-CFc/2.)-CFc/2.)
                if(F(i,k)>0.5) then
                    H1=1.0
                else
                    H1=0.0
                end if 
                if(F(i,k-1)>0.5) then
                    H2=1.0
                else
                    H2=0.
                end if 
                delta=(H1-H2)/dz

                Cur=(W1*Curvc(i,k)+W2*Curvc(i,k-1))/(W1+W2+1.e-20_dp)
                if((2.0*F(i,k)-1.0)*(2.0*F(i,k-1)-1.0)<0.) then
                    Zcsf(i,k)=epiron*delta*Cur
                else
                    Zcsf(i,k)=0.
                end if 
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
    real(kind=dp)::nrx,nrz,sequence1(-1:1,-3:3),sequence2(-3:3,-1:1)
    integer::ii,kk,sgn
    real(kind=dp)::Hr(-1:1),Hz(-1:1)
    real(kind=dp)::Hrzz,Hrz,HZrr,HZr
    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    do i=RS,RE 
        do k=ZS,ZE 
            if(0.<F(i,k).and.F(i,k)<1.) then
                nrx=(F(i+1,k+1)+F(i+1,k-1)-F(i-1,k+1)-F(i-1,k-1)+2.0*(F(i+1,k)-F(i-1,k)))/8./dr
                nrz=(F(i-1,k+1)+F(i+1,k+1)-F(i-1,k-1)-F(i+1,k-1)+2.0*(F(i,k+1)-F(i,k-1)))/8./dz
                if(abs(nrx)<abs(nrz)) then 
                    do ii=-1,1
                        do kk=-3,3
                            sequence1(ii,kk)=F(i+ii,k+kk)
                        end do 
                    end do 
                    call get_Hz(sequence1,nrz,Hz)
                    HZrr=(Hz(1)-2.0*Hz(0)+Hz(-1))/dr**2*dz
                    Hzr=(Hz(1)-Hz(-1))/dr/2.0*dz 
                    Curvc(i,k)=-Hzrr/((sqrt(1.0+Hzr**2))**3+1.0e-20)-Hzr/RC(i)/(sqrt(1.0+Hzr**2)+1.0e-20)
                    if(sgn(nrz*Hzrr)<0.) then
                        Curvc(i,k)=-abs(Curvc(i,k))
                    else
                        Curvc(i,k)=abs(Curvc(i,k))
                    end if 
                else
                    do kk=-1,1
                        do ii=-3,3
                            sequence2(ii,kk)=F(i+ii,k+kk)
                        end do 
                    end do 
                    call get_Hr(sequence2,nrx,Hr)
                    Hrzz=(Hr(1)-2.0*Hr(0)+Hr(-1))/dz**2*dr
                    Hrz=(Hr(1)-Hr(-1))/2.0/dz*dr
                    Curvc(i,k)=-Hrzz/((sqrt(1.0+Hrz**2))**3+1.0e-20)+1.0/RC(i)/(sqrt(1.0+Hrz**2)+1.0e-20)
                    if(sgn(Hrzz*nrx)<0.) then
                        Curvc(i,k)=-abs(Curvc(i,k))
                    else
                        Curvc(i,k)=abs(Curvc(i,k))
                    end if 
                end if 
            end if 
        end do 
    end do

    return

end subroutine curvature_compute

subroutine get_Hz(seq,norin,H)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::seq(-1:1,-3:3)
    real(kind=dp),intent(out)::H(-1:1),norin 
    integer::sgn

    
    if(-sgn(norin)>0.) then
        H(0)=seq(0,-3)+seq(0,-2)+seq(0,-1)+seq(0,0)+&
    seq(0,1)+seq(0,2)+seq(0,3)
        H(-1)=seq(-1,-3)+seq(-1,-2)+seq(-1,-1)+seq(-1,0)+&
    seq(-1,1)+seq(-1,2)+seq(-1,3)
        H(1)=seq(1,-3)+seq(1,-2)+seq(1,-1)+seq(1,0)+&
    seq(1,1)+seq(1,2)+seq(1,3)
    else
        H(0)=7.0-(seq(0,-3)+seq(0,-2)+seq(0,-1)+seq(0,0)+&
    seq(0,1)+seq(0,2)+seq(0,3))
        H(-1)=7.0-(seq(-1,-3)+seq(-1,-2)+seq(-1,-1)+seq(-1,0)+&
    seq(-1,1)+seq(-1,2)+seq(-1,3))
        H(1)=7.0-(seq(1,-3)+seq(1,-2)+seq(1,-1)+seq(1,0)+&
    seq(1,1)+seq(1,2)+seq(1,3))
    end if 

    return

end subroutine get_Hz

subroutine get_Hr(seq,norin,H)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::seq(-3:3,-1:1)
    real(kind=dp),intent(out)::H(-1:1),norin
    integer::sgn

    if(-sgn(norin)>0.) then
        H(0)=seq(-3,0)+seq(-2,0)+seq(-1,0)+seq(0,0)+&
    seq(1,0)+seq(2,0)+seq(3,0)
        H(1)=seq(-3,-1)+seq(-2,-1)+seq(-1,-1)+seq(0,-1)+&
    seq(1,-1)+seq(2,-1)+seq(3,-1)
        H(-1)=seq(-3,1)+seq(-2,1)+seq(-1,1)+seq(0,1)+&
    seq(1,1)+seq(2,1)+seq(3,1)
    else 
        H(0)=7.0-(seq(-3,0)+seq(-2,0)+seq(-1,0)+seq(0,0)+&
    seq(1,0)+seq(2,0)+seq(3,0))
        H(1)=7.0-(seq(-3,-1)+seq(-2,-1)+seq(-1,-1)+seq(0,-1)+&
    seq(1,-1)+seq(2,-1)+seq(3,-1))
        H(-1)=7.0-(seq(-3,1)+seq(-2,1)+seq(-1,1)+seq(0,1)+&
    seq(1,1)+seq(2,1)+seq(3,1))
    end if 

    return

end subroutine get_Hr



integer function  sgn(numbera)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp)::numbera
    if(numbera>0._dp) then
        sgn=1
    else if(numbera<0._dp) then
        sgn=-1
    else
        sgn=0
    end if
end function sgn


