subroutine get_height_function_Hz
	use par 
    implicit none
    integer::iter,ii
    real(kind=dp)::topdown,topup,bottomdown,bottomup
    real(kind=dp)::nz,sign
    real::sgn 
    integer::bas(3)
    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1


   	do i=-1,1
   		bas(i)=0
   	end do 

   	!................................
    do i=RS,RE
            do k=ZS,ZE
            	do ii=-1,1
            		if(0.<F(i+ii,k).and.F(i+ii,k)<1.) then!in the transaint area
            			topdown=F(i+ii,k)
                		bottomup=F(i+ii,k)
                		topup=F(i+ii,k+1)
                		bottomdown=F(i+ii,k-1)
                		Hz(i,k,ii)=F(i+ii,k)*RC(i+ii)*dr*dz  
                		iter=0
                		!start sum from central to base and final 
                		do while(topup/=topdown.or.bottomdown/=bottomup)!one situation exist continue compute
                			if(topup/=topdown) then
                				Hz(i,k,ii)=Hz(i,k,ii)+topup*RC(i+ii)*dr*dz
                			end if 
                			if(bottomdown/=bottomup) then 
                				Hz(i,k,ii)=Hz(i,k,ii)+bottomdown*RC(i+ii)*dr*dz
                			end if 
                     		!valid the base of every height
                     		if(topup>topdown) then  !once base beyond=>up=down
                     			bas(ii)=bas(ii)+1
                     		else if(bottomdown>bottomup) then
                     			bas(ii)=bas(ii)+1
                     		end if 

                			iter=iter+1

                			topup=F(i+ii,k+1+iter)
                			topdown=F(i+ii,k+iter)
                			bottomup=F(i+ii,k-iter)
                			bottomdown=F(i+ii,k-1-iter)

                		end do
                	end if 
                end do !end of ii
                !start compute the relative height
                nz=(F(i,k+1)-F(i,k-1))/dz*0.5_dp
                sign=sgn(nz*(bas(-1)-bas(0)))
                if(abs(bas(-1)-bas(0))==1) then
                    Hz(i,k,-1)=Hz(i,k,-1)-sign*RC(i-1)*dr*dz
                else if(abs(bas(-1)-bas(0))>1) then
                    do l=1,abs(bas(-1)-bas(0))
                	   Hz(i,k,-1)=Hz(i,k,-1)-sign*RC(i-1)*dr*dz
                    end do
                end if  

                sign=sgn(nz*(bas(1)-bas(0)))
                if(abs(bas(1)-bas(0))==1) then
                    Hz(i,k,1)=Hz(i,k,1)-sign*RC(i+1)*dr*dz
                else if(abs(bas(1)-bas(0))>1) then
                    do l=1,abs(bas(1)-bas(0))
                	   Hz(i,k,1)=Hz(i,k,1)-sign*RC(i+1)*dr*dz
                    end do 
                end if 
                !end of correction of relative height
            end do 
    end do 

    return 
end subroutine get_height_function_Hz



subroutine get_height_function_Hr
    use par 
    implicit none
    integer::iter,jj,signal
    real(kind=dp)::righthead,righttail,lefthead,lefttail
    real(kind=dp)::nx,sign
    real::sgn 
    integer::bas(3)
    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1


    do i=-1,1
        bas(i)=0
    end do 

    !................................
    do i=RS,RE
            do k=ZS,ZE
                do jj=-1,1
                    if(0.<F(i,k+jj).and.F(i,k+jj)<1.) then  !in the transiant regine
                        righttail=F(i,k+jj)
                        lefttail=F(i,k+jj)
                        righthead=F(i+1,k+jj)
                        lefthead=F(i-1,k+jj)
                        Hr(i,k,jj)=F(i,k+jj)*RC(i)*dr*dz
                        iter=0
                        do while(righthead/=righttail.and.lefthead/=lefttail)
                            if(righthead/=righttail) then
                                Hr(i,k,jj)=Hr(i,k,jj)+righthead*RC(i+1+iter)*dr*dz
                            end if 
                            if(lefthead/=lefttail) then 
                                Hr(i,k,jj)=Hr(i,k,jj)+lefthead*RC(i-1-iter)*dr*dz
                            end if 

                            if(rihthead>righttail) then  !once base beyond=>head=tail
                                bas(ii)=bas(ii)+1
                                signal=1
                            else if(lefthead>lefttail) then
                                bas(ii)=bas(ii)+1
                                signal=-1
                            end if 

                            iter=iter+1

                            righttail=F(i+iter,k+jj)
                            lefttail=F(i-iter,k+jj)
                            righthead=F(i+1+iter,k+jj)
                            lefthead=F(i-1-iter,k+jj)

                        end do 
                    end if 
                end do !end of jj
                !correction the relative height 
                nx=(F(i,k+1)-F(i,k-1))/dz*0.5_dp
                sign=sgn(nx*(bas(-1)-bas(0)))
                if(abs(bas(-1)-bas(0))==1) then
                    Hr(i,k,-1)=Hr(i,k,-1)-sign*RC(i+signal*(bas(0)+sgn(bas(-1)-bas(0))))*dr*dz
                else if(abs(bas(-1)-bas(0))>1) then
                    do l=1,abs(bas(-1)-bas(0))
                       Hr(i,k,-1)=Hr(i,k,-1)-sign*RC(i+signal*(bas(0)+w*sgn(bas(-1)-bas(0))))*dr*dz
                    end do
                end if  

                sign=sgn(nx*(bas(1)-bas(0)))
                if(abs(bas(1)-bas(0))==1) then
                    Hr(i,k,1)=Hr(i,k,1)-sign*RC(i+signal*(bas(0)+sgn(bas(1)-bas(0))))*dr*dz
                else if(abs(bas(1)-bas(0))>1) then
                    do l=1,abs(bas(1)-bas(0))
                       Hr(i,k,1)=Hr(i,k,1)-sign*RC(i+signal*(bas(0)+w*sgn(bas(1)-bas(0))))*dr*dz
                    end do 
                end if 
            end do 
    end do 

    return 
end subroutine get_height_function_Hr


subroutine compute_curvature
    use par 
    implicit none 

    real(kind=dp)::nx,nz,HZrr,Hzr,Hrzz,Hrz 
 
    call get_height_function_Hz
    call get_height_function_Hr

    RS=1
    ZS=1 
    RE=NR-1
    ZE=NZ-1

    !compute curve
    do i=RS,RE 
        do k=ZS,ZE 
            !reduce the compute complex
            if(F(i,k)/=F(i+1,k).or.F(i,k)/=F(i-1,k).or.F(i,k)/=F(i,k+1).or.F(i,k)/=F(i,k-1)) then 
                nx=(F(i+1,k)-F(i-1,k))/dr*0.5_dp
                nz=(F(i,k+1)-F(ik-1))/dz*0.5_dp
                if(abs(nx)>abs(nz)) then
                    HZrr=(Hz(i,k,1)-2.0_dp*Hz(i,k,0)+Hz(i,k,-1))/dr**2
                    Hzr=(Hz(i,k,1)-Hz(i,k,-1))/dr*0.5_dp
                    Curv(i,k)=-(Hz(i,k,0)**3*Hzrr-Hz(i,k,0)**2*Hzr-Hz(i,k,0)**2)
                    Curv(i,k)=Curv(i,k)/(sqrt(HZ(i,k,0)**2+Hzr**2*Hz(i,k,0)))***3
                else
                    Hrzz=(Hr(i,k,1)-2.0_dp*Hr(i,k,0)+Hr(i,k,-1))/dz**2
                    Hrz=(Hr(i,k,1)-Hr(i,k,-1))/dz*0.5_dp
                    Curv(i,k)=-(Hr(i,k,0)**3*Hrzz-Hr(i,k,0)**2*Hrz-Hr(i,k,0)**2)
                    Curv(i,k)=Curv(i,k)/(sqrt(Hr(i,k,0)**2+Hrz**2*Hr(i,k,0)))**3
                end if 
            else
                Curv(i,k)=0.
            end if 
        end do
    end do

    return 

end subroutine compute_curvature


subroutine compute_csf 
    use par 
    implicit none 

    real(kind=dp)::w1,W2,Cur,delta 

    call compute_curvature

    RS=2
    RE=NR-1
    ZS=1
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                W1=F(i,k)*(1._dp-F(i,k))
                W2=F(i-1,k)*(1._dp-F(i-1,k))

                delta=(F(i,k)-F(i-1,k))/dr

                Cur=(W1*Curv(i,k)+W2*Curv(i-1,k))/(W1+W2+1.e-20_dp)

                 !Rcsf(i,k)=epiron*delta*Cur/(Wdic-Adic)
                Rcsf(i,k)=epiron*delta*Cur
            end do
    end do
    
    
    RS=1
    RE=NR-1
    
    ZS=2
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                W1=F(i,k)*(1._dp-F(i,k))
                W2=F(i,k-1)*(1._dp-F(i,k-1))

                delta=(F(i,k)-F(i,k-1))/dz

                Cur=(W1*Curv(i,k)+W2*Curv(i,k-1))/(W1+W2+1.e-20_dp)

                 !Zcsf(i,k)=epiron*delta*Cur/(Wdic-Adic)
                Zcsf(i,k)=epiron*delta*Cur
            end do
    end do

    return 
end subroutine compute_csf 



