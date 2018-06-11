subroutine iso_vof
    use par 
    implicit none
    
    real(kind=dp)::Vp,pointvof(4),pointloc(4,2)
    real(kind=dp)::Temp1,Temp2,Temptp
    real(kind=dp)::isovalue,cutpos(2,2),x0(2),U0(2),n0(2)
    real(kind=dp)::U0s,facepos(2,2),Tp0(2),Tp(2),FaceAreaInteger
    integer::cell_state,labelt(2),it,iter,itmax
    integer::bound,boundsgn
    integer::cut_test
   
    itmax=10
    !a special condition is that all the pointsvof is equal
    !how to deal with this condition 

    do i=1,NR 
      do k=1,NZ 
        if(i>=2) then
          Fnode(i,k)=(RC(i)*Fo(i,k)+RC(i-1)*Fo(i-1,k)+RC(i)*Fo(i,k-1)+RC(i-1)*Fo(i-1,k-1))
          Fnode(i,k)=Fnode(i,k)/(RC(i)+RC(i-1))/2.0
        else
          Fnode(i,k)=(RC(i)*Fo(i,k)+RC(i)*Fo(i-1,k)+RC(i)*Fo(i,k-1)+RC(i)*Fo(i-1,k-1))
          Fnode(i,k)=Fnode(i,k)/(RC(i)+RC(i))/2.0
        end if 
      end do 
    end do 

    
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1

    !get the fluxs
    !initial the vofvolumes on faces
    !with upwind values
    do i=1,NR 
        do k=1,NZ-1
            fluxr(i,k)=U(i,k)*RP(i)*dz
            if(fluxr(i,k)>0.) then
                vof_face_volr(i,k)=Fo(i-1,k)*fluxr(i,k)*dt
            else
                vof_face_volr(i,k)=Fo(i,k)*fluxr(i,k)*dt
            end if 
        end do 
    end do 

    do k=1,NZ
        do i=1,NR-1
            fluxz(i,k)=W(i,k)*RC(i)*dr 
            if(fluxz(i,k)>0.) then
                vof_face_volz(i,k)=Fo(i,k-1)*fluxz(i,k)*dt
            else
                vof_face_volz(i,k)=Fo(i,k)*fluxz(i,k)*dt
            end if 
        end do 
    end do 

    !step 2 update dvf for interface cells
    do i=RS,RE 
        do k=ZS,ZE
            if(cell_state(Fo(i,k),10*mass)==0) then
                labelt(1)=i
                labelt(2)=k
                !basic information will be used in 
                !next step
                call nodesvof(labelt,Pointvof)
            
                call vertexloc(labelt,Pointloc)

                !call cutpolygeo(F(i,k),Pointvof,Pointloc,isovalue)
                cut_test=0
      
                call vofcutcell(cut_test,Fo(i,k),Pointvof,Pointloc,isovalue)
            
                call isoinfo(isovalue,Pointvof,Pointloc,cutpos,x0,n0)
                call velinterp(labelt,x0,U0)
                !moving velocity(u0*n0)
                U0s=U0(1)*n0(1)+U0(2)*n0(2)

                !two condition is exit
                !if(U0s>0.) then
                  !Ufl is the component of U0s on x,z direction
                  !if U0,x0 in the same direction,Ufl=(U0*n0)n0
                  !if U0,x0 in the opposite direction Ufl=(U0*-n0)(-n0)
                  !......................w  face
                  !update when it is a downwind face 
                  if(fluxr(i,k)<0.) then
                    Facepos(1,1)=RP(i)
                    Facepos(1,2)=ZP(k)
                    Facepos(2,1)=RP(i)
                    Facepos(2,2)=ZP(k+1)
                    !(P-x0)*n0/U0s
                    do it=1,2
                        Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                    end do
                    call order(Tp0,2,Tp) 
                    if(Tp0(1)>Tp0(2)) then
                      Temp1=Facepos(1,1)
                      Temp2=Facepos(1,2)
                      Facepos(1,1)=Facepos(2,1)
                      Facepos(1,2)=Facepos(2,2)
                      Facepos(2,1)=Temp1
                      Facepos(2,2)=Temp2
                      Temptp=Tp0(1)
                      Tp0(1)=TP0(2)
                      Tp0(2)=Temptp
                    end if 
                    call timeintegerArea(.true.,Tp,dt,facepos,FaceAreaInteger)
                    if(U0s<0) then
                      vof_face_volr(i,k)=fluxr(i,k)*dt-vof_face_volr(i,k)
                    end if 
                    if(abs(U(i,k)*FaceAreaInteger)<abs(vof_face_volr(i,k))) then
                        vof_face_volr(i,k)=U(i,k)*FaceAreaInteger
                        if(abs(U(i,k)*FaceAreaInteger-vof_face_volr(i,k))>1.0e-8) then
                          write(*,*)'compute false--west'
                        end if 
                    end if 
                    !the light fluid is moving to the heavy fuild
                    if(U0s<0) then
                      vof_face_volr(i,k)=fluxr(i,k)*dt-vof_face_volr(i,k)
                    end if 
                  end if 
                  !...........................e face 
                  if(fluxr(i+1,k)>0.) then
                    Facepos(1,1)=RP(i+1)
                    Facepos(1,2)=ZP(k)
                    Facepos(2,1)=RP(i+1)
                    Facepos(2,2)=ZP(k+1)
                    !(P-x0)*n0/U0s
                    do it=1,2
                        Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                    end do 

                    call order(Tp0,2,Tp)
                    if(Tp0(1)>Tp0(2)) then
                      Temp1=Facepos(1,1)
                      Temp2=Facepos(1,2)
                      Facepos(1,1)=Facepos(2,1)
                      Facepos(1,2)=Facepos(2,2)
                      Facepos(2,1)=Temp1
                      Facepos(2,2)=Temp2
                      Temptp=Tp0(1)
                      Tp0(1)=TP0(2)
                      Tp0(2)=Temptp
                    end if
                    call timeintegerArea(.true.,Tp,dt,facepos,FaceAreaInteger)
                
                    if(abs(U(i+1,k)*FaceAreaInteger)<abs(vof_face_volr(i+1,k))) then
                        vof_face_volr(i+1,k)=U(i+1,k)*FaceAreaInteger
                        if(abs(U(i+1,k)*FaceAreaInteger-vof_face_volr(i+1,k))>1.0e-8) then
                          write(*,*)'compute false--east'
                        end if 
                    end if 
                     !the light fluid is moving to the heavy fuild
                    if(U0s<0) then
                      vof_face_volr(i+1,k)=fluxr(i+1,k)*dt-vof_face_volr(i+1,k)
                    end if

                  end if 
                  !............................bottom face
                  if(fluxz(i,k)<0.) then
                    Facepos(1,1)=RP(i)
                    Facepos(1,2)=ZP(k)
                    Facepos(2,1)=RP(i+1)
                    Facepos(2,2)=ZP(k)
                  !(P-x0)*n0/U0s
                    do it=1,2
                        Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                    end do 
                    call order(Tp0,2,Tp)
                    if(Tp0(1)>Tp0(2)) then
                      Temp1=Facepos(1,1)
                      Temp2=Facepos(1,2)
                      Facepos(1,1)=Facepos(2,1)
                      Facepos(1,2)=Facepos(2,2)
                      Facepos(2,1)=Temp1
                      Facepos(2,2)=Temp2
                      Temptp=Tp0(1)
                      Tp0(1)=TP0(2)
                      Tp0(2)=Temptp
                    end if
                    call timeintegerArea(.false.,Tp,dt,facepos,FaceAreaInteger)

                    if(abs(W(i,k)*FaceAreaInteger)<abs(vof_face_volz(i,k))) then
                        vof_face_volz(i,k)=W(i,k)*FaceAreaInteger
                        if(abs(W(i,k)*FaceAreaInteger-vof_face_volz(i,k))>1.0e-8) then
                          write(*,*)'compute false--south'
                        end if 
                    end if 

                    !the light fluid is moving to the heavy fuild
                    if(U0s<0) then
                      vof_face_volz(i,k)=fluxz(i,k)*dt-vof_face_volz(i,k)
                    end if 

                  end if 
                  !..............................top face
                  if(fluxz(i,k+1)>0.) then
                    Facepos(1,1)=RP(i)
                    Facepos(1,2)=ZP(k+1)
                    Facepos(2,1)=RP(i+1)
                    Facepos(2,2)=ZP(k+1)
                    !(P-x0)*n0/U0s
                    do it=1,2
                        Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                    end do 
                    call order(Tp0,2,Tp)
                    if(Tp0(1)>Tp0(2)) then
                      Temp1=Facepos(1,1)
                      Temp2=Facepos(1,2)
                      Facepos(1,1)=Facepos(2,1)
                      Facepos(1,2)=Facepos(2,2)
                      Facepos(2,1)=Temp1
                      Facepos(2,2)=Temp2
                      Temptp=Tp0(1)
                      Tp0(1)=TP0(2)
                      Tp0(2)=Temptp
                    end if
                    call timeintegerArea(.false.,Tp,dt,facepos,FaceAreaInteger)

                    if(abs(W(i,k+1)*FaceAreaInteger)<abs(vof_face_volz(i,k+1))) then
                        vof_face_volz(i,k+1)=W(i,k+1)*FaceAreaInteger
                        if(abs(W(i,k+1)*FaceAreaInteger-vof_face_volz(i,k+1))>1.0e-8) then
                          write(*,*)'compute false--north'
                        end if 
                    end if 

                    !the light fluid is moving to the heavy fuild
                    if(U0s<0) then
                      vof_face_volz(i,k+1)=fluxz(i,k+1)*dt-vof_face_volz(i,k+1)
                    end if 

                  end if 
                  !end the four condition

            end if 
        end do 
    end do 

    !initial the bound state
    do i=RS,RE 
      do k=ZS,ZE 
        check_bound(i,k)=.false.
      end do 
    end do 
    !at first set the bound as true
    bound=1
    !correct the volume fractions
    !and correct the bound state 
    do i=RS,RE 
      do k=ZS,ZE 
              Vp=RC(i)*dr*dz
              F(i,k)=Fo(i,k)-(-vof_face_volr(i,k)+vof_face_volr(i+1,k)-vof_face_volz(i,k)+vof_face_volz(i,k+1))/Vp

!             !1.update
!             !2.correct itself 
!             !3.call it's neighbor to complete 1 and 2
!             call correct_recusion(i,k)
      end do 
    end do 
    !call boundary F(i,k)
    call boundary_f

    do i=RS,RE 
      do k=ZS,ZE 
        if(F(i,k)<0..or.F(i,k)>1.) then
        

          boundsgn=0
          check_bound(i,k)=.false.
        else
          boundsgn=1
        end if 
        bound=bound*boundsgn
      end do 
    end do 

    !step3 correct volumes for the unbound cells
    do while(bound==0.and.iter<=itmax)
        bound=1

        do i=RS,RE
           do k=ZS,ZE
              call correct_recusion(i,k)
            end do 
        end do 
        call boundary_f

        do i=RS,RE 
          do k=ZS,ZE 
            if(F(i,k)<0..or.F(i,k)>1.) then
              boundsgn=0
            else
              boundsgn=1
            end if 
            bound=bound*boundsgn
          end do 
        end do 
        !end of correct 

        iter=iter+1
    end do !end of correct


   !finial step:set tolerance for compute 
    do i=RS,RE 
      do k=ZS,ZE 
          if(F(i,k)<0..and.F(i,k)>-mass) then
                F(i,k)=0.
          else if(F(i,k)>1..and.F(i,k)-1.<mass) then
                F(i,k)=1.
          end if 
          if(F(i,k)<0..or.F(i,k)>1.) then
            write(*,*)'false exit'
            write(*,*)F(i,k),i,k
          end if 
      end do 
    end do 
    call boundary_f

    return
end subroutine iso_vof

subroutine timeintegerArea(isrdirectionface,Tpin,deltat,faceloc,IntegerArea)
	!integer area
  !dvf=S(tao)*dtao
  !Vf=integer(dvf)
  !iszdirection dvf=U0(2)*(tao-tp1)*dtao
  !iszdirection dvf=0.5*((u0(1)*(tao-tp1)+xp1)**2-xp1**2)*dtao
  implicit none
	integer,parameter::dp=selected_real_kind(p=15)
	logical,intent(in)::isrdirectionface
	real(kind=dp),intent(in)::Tpin(2),deltat
	real(kind=dp),intent(in)::faceloc(2,2)
	real(kind=dp),intent(out)::IntegerArea
	real(kind=dp)::IntegerArea1,IntegerArea2,isoU0(2)
	if(isrdirectionface) then !loc(1,1)=loc(2,1) Tpin(1)<Tpin(2)
		!r1*(z2-z1) nf is in r direction
    !defalut setting:Tpin(1)>Tpin(2),
    isoU0(2)=(faceloc(2,2)-faceloc(1,2))/(Tpin(2)-Tpin(1))
		if(Tpin(1)>deltat) then
			IntegerArea=0.
		else if(Tpin(2)<0.) then
			IntegerArea=abs(faceloc(1,1)*(faceloc(2,2)-faceloc(1,2))*deltat)
		else if(Tpin(1)<0..and.Tpin(2)>0.) then
      !tp2
      !deltat
      !cut==t=0
      !tp1
			if(Tpin(2)>=deltat) then
				IntegerArea1=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(0.-Tpin(1))**2)
				IntegerArea2=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(deltat-Tpin(1))**2)
				IntegerArea=IntegerArea2-IntegerArea1
			else if(Tpin(2)<deltat) then
      !deltat
      !tp2
      !cut==t=0
      !tp1
				IntegerArea1=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(0.-Tpin(1))**2)
				IntegerArea2=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(Tpin(2)-Tpin(1))**2)
				IntegerArea=IntegerArea2-IntegerArea1
				IntegerArea=IntegerArea+abs(faceloc(1,1)*(faceloc(2,2)-faceloc(1,2))*(deltat-Tpin(2)))
			end if 
    else if(Tpin(1)>0..and.Tpin(2)>0.) then
      !tp2
      !deltat
      !tp1
      !cut=t=0.
      if(Tpin(2)>=deltat) then
        IntegerArea1=0.
        IntegerArea2=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(deltat-Tpin(1))**2)
        IntegerArea=IntegerArea2-IntegerArea1
      !deltat
      !tp2
      !tp1
      !cut=t=0.
      else if(Tpin(2)<deltat) then
        IntegerArea1=0.
        IntegerArea2=abs(0.5_dp*isoU0(2)*faceloc(1,1)*(Tpin(2)-Tpin(1))**2)
        IntegerArea=IntegerArea2-IntegerArea1
        IntegerArea=IntegerArea+abs(faceloc(1,1)*(faceloc(2,2)-faceloc(1,2))*(deltat-Tpin(2)))
      end if 
		end if 
	else
		!zdirectionface nf is in z direction
		!Sface=(r1+r1)/2*(r2-r1)
    isoU0(1)=(faceloc(2,1)-faceloc(1,1))/(Tpin(2)-Tpin(1))
		if(Tpin(1)>deltat) then
			IntegerArea=0.
		else if(Tpin(2)<0.) then
			IntegerArea=abs((faceloc(1,1)+faceloc(2,1))*0.5_dp*(faceloc(2,1)-faceloc(1,1))*deltat)
		else if(Tpin(1)<0..and.Tpin(2)>0.) then
			if(Tpin(2)>deltat) then
        IntegerArea1=abs(1./6.*isoU0(1)**2*(0.-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(0.-Tpin(1))**2)
        IntegerArea2=abs(1./6.*isoU0(1)**2*(deltat-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(deltat-Tpin(1))**2)
				IntegerArea=IntegerArea2-IntegerArea1
			else if(Tpin(2)<deltat) then
				IntegerArea1=abs(1./6.*isoU0(1)**2*(0.-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(0.-Tpin(1))**2)
				IntegerArea2=abs(1./6.*isoU0(1)**2*(Tpin(2)-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(Tpin(2)-Tpin(1))**2)
				IntegerArea=IntegerArea2-IntegerArea1
				IntegerArea=IntegerArea+abs((faceloc(1,1)+faceloc(2,1))*0.5_dp*(faceloc(2,1)-faceloc(1,1))*(deltat-Tpin(2)))
			end if 
    else if(Tpin(1)>0..and.Tpin(2)>0.) then 
      if(Tpin(2)>=deltat) then
        !0-Tpin(1) area=0
        !Tpin(1)-deltat
        IntegerArea1=0.
        IntegerArea2=abs(1./6.*isoU0(1)**2*(deltat-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(deltat-Tpin(1))**2)
        IntegerArea=IntegerArea2-IntegerArea1
      else if(Tpin(2)<deltat) then
        IntegerArea1=0.
        IntegerArea2=abs(1./6.*isoU0(1)**2*(Tpin(2)-Tpin(1))**3+0.5*faceloc(1,1)*isoU0(1)*(Tpin(2)-Tpin(1))**2)
        IntegerArea=IntegerArea2-IntegerArea1
        IntegerArea=IntegerArea+abs((faceloc(1,1)+faceloc(2,1))*0.5_dp*(faceloc(2,1)-faceloc(1,1))*(deltat-Tpin(2)))
      end if 
		end if 
	end if 
	return
end subroutine timeintegerArea

!no use of par 
!subroutine face_volume_correct has no conflict with par 

recursive subroutine correct_recusion(la,lb)
    !almost every varible are global varibales
    use par 
    implicit none 
    integer::it,correct_test
    integer,intent(in)::la,lb 
    real(kind=dp)::Phi(4),face_vlome(4),face_vlomeout(4),volumesuper
    real(kind=dp)::Vp 
    !update itself
    Vp=dr*dz*Rc(la)
    F(la,lb)=Fo(la,lb)-(-vof_face_volr(la,lb)+vof_face_volr(la+1,lb)-vof_face_volz(la,lb)+vof_face_volz(la,lb+1))/Vp
    !basic procession
    if(0<=F(la,lb).and.F(la,lb)<=1.) then
      check_bound(la,lb)=.true.
    else

        !correct itself 
        correct_test=0

        if(F(la,lb)>1.) then
            volumesuper=(F(la,lb)-1.)*Vp
            do it=1,2
            !convert to positive values
                Phi(it)=(-1)**it*fluxr(la+it-1,lb)
                Phi(it+2)=(-1)**it*fluxz(la,lb+it-1)
                face_vlome(it)=(-1)**it*vof_face_volr(la+it-1,lb)
                face_vlome(it+2)=(-1)**it*vof_face_volz(la,lb+it-1)
            end do 
         
            call face_volume_correct(correct_test,Phi,face_vlome,volumesuper,dt,face_vlomeout)
            do it=1,2
              !convert to original values
                vof_face_volr(la+it-1,lb)=(-1)**it*face_vlomeout(it)
                vof_face_volz(la,lb+it-1)=(-1)**it*face_vlomeout(it+2)
            end do 
        else if(F(la,lb)<0.) then
            !do some convertion,apla=>a-alpha,dVf=Phi*dt-dVf
            volumesuper=(-F(la,lb))*Vp
            do it=1,2
            !convert to positive values
                Phi(it)=(-1)**it*fluxr(la+it-1,lb)
                Phi(it+2)=(-1)**it*fluxz(la,lb+it-1)
                face_vlome(it)=(-1)**it*(fluxr(la+it-1,lb)*dt-vof_face_volr(la+it-1,lb))
                face_vlome(it+2)=(-1)**it*(fluxz(la,lb+it-1)*dt-vof_face_volz(la,lb+it-1))
            end do 
     
            call face_volume_correct(correct_test,Phi,face_vlome,volumesuper,dt,face_vlomeout)
            do it=1,2
              !convert to original values
                vof_face_volr(la+it-1,lb)=fluxr(la+it-1,lb)*dt-(-1)**it*face_vlomeout(it)
                vof_face_volz(la,lb+it-1)=fluxz(la,lb+it-1)*dt-(-1)**it*face_vlomeout(it+2)
            end do 
        end if 

        F(la,lb)=Fo(la,lb)-(-vof_face_volr(la,lb)+vof_face_volr(la+1,lb)-vof_face_volz(la,lb)+vof_face_volz(la,lb+1))/Vp
        check_bound(la,lb)=.true.




        !for it's downwind neighbor,call for a recusion
        !do not need to set_boundary condition,for there is no 
        !downwind cells near the boundary_face
        do it=1,2
          if(face_vlomeout(it)>0.) then
            !boundary condition 
            if(la+(-1)**it>=1.and.la+(-1)**it<=NR-1) then
                call correct_recusion(la+(-1)**it,lb)
                
            end if 
          end if 
          if(face_vlomeout(it+2)>0.) then
            !boundary_condition
            if(lb+(-1)**it>=1.and.lb+(-1)**it<=NZ-1) then
                call correct_recusion(la,lb+(-1)**it)
            end if
          end if 
        end do 
        !check_over all the neighbors
    end if 
    return
end subroutine correct_recusion



subroutine face_volume_correct(testnumber,Phiin,F1volume,volumeplus,deltat,F1volumeout)
    !correct the volume pass the face during [t,t+\delta t]
    implicit none
    integer,intent(In)::testnumber
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp)::Phiin(4),F1volume(4)
    logical::isdownwind(4)
    real(kind=dp)::volumeplus,deltat,testvol
    real(kind=dp),intent(out)::F1volumeout(4)
    real(kind=dp)::total,totaltemp,deliveredvolume,a
    integer::deliversize,it,i 
    logical::capacity(4)

    !decide the values of ifisdownwind
    do i=1,4
        if(Phiin(i)>0.) then
            isdownwind(i)=.true.
        else
            isdownwind(i)=.false.
        end if 
    end do 
    !start
    total=0.
    deliversize=0
    do i=1,4
        if(isdownwind(i).and.F1volume(i)<Phiin(i)*deltat) then
            capacity(i)=.true.
            total=total+Phiin(i)
            deliversize=deliversize+1
        else
            capacity(i)=.false.
        end if 
    end do 

    !convert at out,volumplus is setting >0.
    !start deliver
    totaltemp=total 
    it=1
    do while(volumeplus>0..and.it<=deliversize) 
            !volumeplus==0. deliver over
            !it==deliversize,all downwindface get its capacity
            deliveredvolume=0.

            !start a new deliver with (vloumeplus,capacity,total)
            do i=1,4
                if(capacity(i)) then
                    !recorrect the capacity ability
                    if(volumeplus*Phiin(i)/total+F1volume(i)>=Phiin(i)*deltat) then 
                        !stitl use the flux as weight
                        !use the totaltemp to save the capicity change
                        totaltemp=totaltemp-Phiin(i)
                        capacity(i)=.false. !capacity reached
                        a=Phiin(i)*deltat-F1volume(i)
                        deliveredvolume=deliveredvolume+a
                        F1volume(i)=Phiin(i)*deltat
                    else
                        a=volumeplus*Phiin(i)/total
                        deliveredvolume=deliveredvolume+a
                        F1volume(i)=F1volume(i)+a 
                    end if 
                end if 
            end do 
            !the rest volumeplus
            volumeplus=volumeplus-deliveredvolume
            !change the total after i iteration from 1 to 4
            total=totaltemp
            it=it+1
  end do 
  !convert over 
  do i=1,4
        F1volumeout(i)=F1volume(i)
  end do 

  return
end subroutine  face_volume_correct

subroutine vofcutcell(test_number,cellvof,vertexvof,vertexloc,f3)
	!expanation: !get the iso cut value to make a(f) = alpha1
  implicit none
  integer,parameter::dp=selected_real_kind(p=15)
  integer,intent(in)::test_number
  real(kind=dp),intent(in)::cellvof,vertexvof(4),vertexloc(4,2)
  real(kind=dp),intent(out)::f3
  real(kind=dp)::res,tol,small,alpha1,alpha
  real(kind=dp)::forder(4),a1,a2,a3,a4,f1,f2,f4
  real(kind=dp)::x0,x1,x2,g0,g1,g2
  real(kind=dp)::a(4),fi(4),M(4,4),C(4)
  integer::itf,itfmax,L1,L2,L3,icy,jcy
  logical::compute_wella1,compute_wella2

  tol=1.e-10
  small=1.e-6
  itfmax=100
  
  alpha1=cellvof 
	!Finding the two vertices inbetween which the isovalue giving alpha1 lies
	!get a cell's nodes values//Fnodes
	!sorted the array


	call order(vertexvof,4,forder)

	f1=forder(1)
	f2=forder(4)

	L1=1
	L2=4

	a1=1. !f1N=0
	a2=0. !f1n=1
  !the guess a1,a2 a1,a2 must be accurate compute
  compute_wella1=.false.
  compute_wella2=.false.
  !assume decrease with the increase of f 
  !make sure that alpha1 in [f1,f2]
	do while(L2-L1>1)
		L3=int(0.5_dp*(L1+L2))
		f3=forder(L3)
    call cutpolygeo(f3,vertexvof,vertexloc,a3)

		if(a3>alpha1) then
      L1=L3
      f1=f3
      a1=a3
      compute_wella1=.true.
		else if(a3<alpha1) then
      L2=L3
      f2=f3 
      a2=a3
      compute_wella2=.true.
		end if 
	end do 
  !if a1,a2 is not compute accurate
  if(.not.(compute_wella1)) then
    call cutpolygeo(f1,vertexvof,vertexloc,a1)
  end if 

  if(.not.(compute_wella2)) then
    call cutpolygeo(f2,vertexvof,vertexloc,a2)
  end if 

  if(abs(f2-f1)<10*small) then
      f3=f1
  else if(abs(a2-a1)<tol) then
      f3=0.5_dp*(f1+f2)
  else
	     !Now we know that a(f) = alpha1 is to be found on the f interval
      ![f1, f2], i.e. alpha1 will be in the interval [a2,a1]
      !Finding coefficients in 3 deg polynomial alpha(f) from 4 solutions
      f3=f1+1./3.*(f2-f1)
      call cutpolygeo(f3,vertexvof,vertexloc,a3)
      f4=f1+2./3.*(f2-f1)
      call cutpolygeo(f4,vertexvof,vertexloc,a4)
      !Building and solving Vandermonde matrix equation
      !M*C=a
   	  a(1)=a2  !f2
   	  a(2)=a4   
   	  a(3)=a3
   	  a(4)=a1  !f1
   	  fi(1)=1.
   	  fi(2)=(f4-f1)/(f2-f1)
   	  fi(3)=(f3-f1)/(f2-f1)
   	  fi(4)=0.

   	  do icy=1,4
   		   do jcy=1,4
   			    M(icy,jcy)=fi(icy)**(4-jcy)
   		   end do 
   	  end do 
      !M*C=a return C
      !(left n*n,right b,n,cof)

   	  call LUsolver(M,a,4,C)

   	  !find root with Newton method
   	  !https://en.wikipedia.org/wiki/Newton%27s_method
   	  !congver from biggest one to the finial results
   	  f3=fi(1) !1 give the maxvalue to avoid the singularity//f3=0,C(3)=0
   	  a3=a(1)
   	  itf=0
   	  res=abs(a3-alpha1)
   	  do while(res>tol.and.itf<=itfmax) 
   		   f3=f3-(C(1)*f3**3 + C(2)*f3**2 + C(3)*f3 + C(4) - alpha1)/(3*C(1)*f3**2 + 2*C(2)*f3 + C(3))
   		   a3=C(1)*f3**3 + C(2)*f3**2 + C(3)*f3 + C(4)
   		   res=abs(a3-alpha1)
   		   itf=itf+1
   	  end do 

   	  !Scaling back to original range
      f3 = f3*(f2 - f1) + f1
      !check results 
      call cutpolygeo(f3,vertexvof,vertexloc,alpha)
      res=abs(alpha-alpha1)

      !if do not meet,use Secant method to found the root
      !https://en.wikipedia.org/wiki/Secant_method
      !for the feaures of Newton method,even do not converge
      !for the Secant,We do not need a clear function
      !make sure the the conveger value is in middle

      if(res>tol) then
    	   !x2=f3 !the low limited
    	   !g2=alpha-alpha1
         !x1=max(1e-3*(f2-f1),100*small)
         !x1 = max(x1,f1)
         !x1 = min(x1,f2)
    	   !call cutpolygeo(x1,vertexvof,vertexloc,alpha)
    	   !g1=alpha-alpha1
         x1 = f1
         g1 = a1 - alpha1
         x2 = f2
         g2 = a2 - alpha1
       
    	   itf=0
    	   do while(res>tol.and.itf<=itfmax.and.g1/=g2)
    		      x0=(x2*g1 - x1*g2)/(g1 - g2)
    		      call cutpolygeo(x0,vertexvof,vertexloc,alpha)
    		      g0=alpha-alpha1
    		      res=abs(g0)
    		      x2 = x1
    		      g2 = g1
              x1 = x0
              g1 = g0

              itf=itf+1
         end do 
         f3=x0 !now assume that we have complete the search 
      end if 

  end if 

  return 
end subroutine vofcutcell 


  
    