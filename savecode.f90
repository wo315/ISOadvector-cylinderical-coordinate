
                    if(abs(U(i,k)*FaceAreaInteger)<abs(vof_face_volr(i,k))) then
                        vof_face_volr(i,k)=U(i,k)*FaceAreaInteger
                        if(abs(U(i,k)*FaceAreaInteger-vof_face_volr(i,k))>1.0e-8) then
                          write(*,*)'compute false--west'
                        end if 
                    end if                 if(i==4.and.k==31) then
                  write(*,*)"Pointloc(1,1),pointloc(1,2)"
                  write(*,*)Pointloc(1,1),pointloc(1,2)
                  write(*,*)"Pointloc(2,1),pointloc(2,2)"
                  write(*,*)Pointloc(2,1),pointloc(2,2)
                  write(*,*)"Pointloc(3,1),pointloc(3,2)"
                  write(*,*)Pointloc(3,1),pointloc(3,2)
                  write(*,*)"Pointloc(4,1),pointloc(4,2)"
                  write(*,*)Pointloc(4,1),pointloc(4,2)
                  write(*,*)"cutpos(1,1),cutpos(1,2)"
                  write(*,*)cutpos(1,1),cutpos(1,2)
                 write(*,*)"cutpos(2,1),cutpos(2,2)"
                  write(*,*)cutpos(2,1),cutpos(2,2)
                end if 
        if(la==40.and.lb==40) then
            write(*,*)'correct befor'
            write(*,*)F(la,lb),FO(la,lb),la,lb
            write(*,*)fluxr(la,lb)*dt,fluxr(la+1,lb)*dt
            write(*,*)vof_face_volr(la,lb),vof_face_volr(la+1,lb)
            write(*,*)fluxz(la,lb)*dt,fluxz(la,lb+1)*dt
            write(*,*)vof_face_volz(la,lb),vof_face_volz(la,lb+1)
        endif
          if(i==11.and.k==37) then
            write(*,*)F(i,k),FO(i,k)
            write(*,*)fluxr(i,k)*dt,fluxr(i+1,k)*dt
            write(*,*)vof_face_volr(i,k),vof_face_volr(i+1,k)
            write(*,*)fluxz(i,k)*dt,fluxz(i,k+1)*dt
            write(*,*)vof_face_volz(i,k),vof_face_volz(i,k+1)
          end if 
        if(i==91.and.k==155) then
            write(*,*)'correct after'
            write(*,*)F(i,k),FO(i,k),i,k
            write(*,*)fluxr(la,lb)*dt,fluxr(la+1,lb)*dt
            write(*,*)vof_face_volr(la,lb),vof_face_volr(la+1,lb)
            write(*,*)fluxz(la,lb)*dt,fluxz(la,lb+1)*dt
            write(*,*)vof_face_volz(la,lb),vof_face_volz(la,lb+1)
        endif

   subroutine vof_update
    use par 
    implicit none
    real(kind=dp),pointer,dimension(:,:)::fluxr,fluxz,vof_face_volr,vof_face_volz
    real(kind=dp)::Vp,pointvof(4),pointloc(4,2)
    real(kind=dp)::isovalue,cutpos(2,2),x0(2),U0(2),n0(2)
    real(kind=dp)::U0s,facepos(2,2),Tp0(2),Tp(2),FaceAreaInteger
    real(kind=dp)::Phi(4),face_vlome(4),face_vlomeout(4),volumesuper
    integer::cell_state,labelt(2),it,iter,itmax
    logical::bound

    allocate(fluxr(1:NR,1:NZ),fluxz(1:NR,1:NZ),vof_face_volr(1:NR,1:NZ),vof_face_volz(1:NR,1:NZ))
   
    itmax=10
    
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
                vof_face_volr(i,k)=F(i-1,k)*fluxr(i,k)*dt
            else
                vof_face_volr(i,k)=F(i,k)*fluxr(i,k)*dt
            end if 
        end do 
    end do 

    do k=1,NZ
        do i=1,NR-1
            fluxz(i,k)=W(i,k)*RC(i)*dr 
            if(fluxz(i,k)>0.) then
                vof_face_volz(i,k)=F(i,k-1)*fluxz(i,k)*dt
            else
                vof_face_volz(i,k)=F(i,k)*fluxz(i,k)*dt
            end if 
        end do 
    end do 
    !step 2 update dvf for interface cells
    do i=RS,RE 
        do k=ZS,ZE
            Vp=dr*dz*Rc(i)
            if(cell_state(F(i,k))==0) then
                labelt(1)=i
                labelt(2)=k
                !basic information will be used in 
                !next step
                call nodesvof(labelt,Pointvof)
                call vertexloc(labelt,Pointloc)
                !call cutpolygeo(F(i,k),Pointvof,Pointloc,isovalue)
                call vofcutcell(labelt,Pointvof,Pointloc,isovalue)
                !call isoinfo(isovalue,Pointvof,Pointloc,cutpos,x0,n0)
                !call velinterp(labelt,x0,U0)
                !moving velocity(u0*n0)
                U0s=U0(1)*n0(1)+U0(2)*n0(2)
                !......................w  face
                Facepos(1,1)=RP(i)
                Facepos(1,2)=ZP(k)
                Facepos(2,1)=RP(i)
                Facepos(2,2)=ZP(k+1)
                !(P-x0)*n0/U0s
                do it=1,2
                    Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                end do 
                call order(Tp0,2,Tp)
                call timeintegerArea(.true.,Tp,dt,facepos,U0,FaceAreaInteger)
                vof_face_volr(i,k)=U(i,k)*FaceAreaInteger
                !...........................e face 
                Facepos(1,1)=RP(i+1)
                Facepos(1,2)=ZP(k)
                Facepos(2,1)=RP(i+1)
                Facepos(2,2)=ZP(k+1)
                !(P-x0)*n0/U0s
                do it=1,2
                    Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                end do 
                call order(Tp0,2,Tp)
                call timeintegerArea(.true.,Tp,dt,facepos,U0,FaceAreaInteger)
                vof_face_volr(i+1,k)=U(i+1,k)*FaceAreaInteger
                !............................bottom face
                Facepos(1,1)=RP(i)
                Facepos(1,2)=ZP(k)
                Facepos(2,1)=RP(i+1)
                Facepos(2,2)=ZP(k)
                !(P-x0)*n0/U0s
                do it=1,2
                    Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                end do 
                call order(Tp0,2,Tp)
                call timeintegerArea(.false.,Tp,dt,facepos,U0,FaceAreaInteger)
                vof_face_volz(i,k)=W(i,k)*FaceAreaInteger
                !..............................top face
                Facepos(1,1)=RP(i)
                Facepos(1,2)=ZP(k+1)
                Facepos(2,1)=RP(i+1)
                Facepos(2,2)=ZP(k+1)
                !(P-x0)*n0/U0s
                do it=1,2
                    Tp0(it)=((facepos(it,1)-x0(1))*n0(1)+(facepos(it,2)-x0(2))*n0(2))/U0s
                end do 
                call order(Tp0,2,Tp)
                call timeintegerArea(.false.,Tp,dt,facepos,U0,FaceAreaInteger)
                vof_face_volz(i,k+1)=W(i,k+1)*FaceAreaInteger
            end if 

            !correct for the boundary conditions
            if(i==RS.or.i==NR) then
                vof_face_volr(i,k)=0.
            end if 
            if(k==ZS.or.k==NZ) then
                vof_face_volz(i,k)=0.
            end if 
            !no all the vof_face has been compute
            !update the F for the next time step
            F(i,k)=F(i,k)-(vof_face_volr(i,k)+vof_face_volr(i+1,k)+vof_face_volz(i,k)+vof_face_volz(i,k+1))/Vp
        end do 
    end do 
    
    !call boundary F(i,k)
    call boundary_f
    !step3 correct volumes for the unbound cells
    do while(bound.eqv..false..and.iter<=itmax)
        do i=RS,RE
           do k=ZS,ZE
                Vp=dr*dz*RC(i)
                if(F(i,k)<0..or.F(i,k)>1.) then
                     do it=1,2
                        Phi(it)=fluxr(i+it-1,k)
                        Phi(it+2)=fluxz(i,k+it-1)
                        face_vlome(it)=vof_face_volr(i+it-1,k)
                        face_vlome(it+2)=vof_face_volz(i,k+it-1)
                     end do 
                     if(F(i,k)<0.) then
                        volumesuper=F(i,k)*Vp
                     else if(F(i,k)>1.) then
                        volumesuper=(F(i,k)-1.)*Vp
                     end if 
                     call face_volume_correct(Phi,face_vlome,volumesuper,dt,face_vlomeout)
                     do it=1,2
                        vof_face_volr(i+it-1,k)=face_vlomeout(it)
                        vof_face_volz(i,k+it-1)=face_vlomeout(it+2)
                     end do 
                end if 
                !reupdate the vof values
                F(i,k)=F(i,k)-(vof_face_volr(i,k)+vof_face_volr(i+1,k)+vof_face_volz(i,k)+vof_face_volz(i,k+1))/Vp
                if(F(i,k)<0..or.F(i,k)>1.) then
                    bound=.false.
                else
                    bound=.true.
                end if 
           end do 
        end do 
        call boundary_f
        iter=iter+1
    end do !end of correct

    deallocate(fluxr,fluxz,vof_face_volr,vof_face_volz)
    return
end subroutine vof_update
    do itout=0,1
        do itin=0,1
            do i=0,NR
                do k=0,NZ
                    smoothvalue1(i,k,0,0)=Nnormalcr(i,k) 
                    smoothvalue2(i,k,0,0)=Nnormalcz(i,k) 
                end do 
            end do 

            !start compute
            do i=RS,RE
                do k=ZS,ZE
                    if(itout==0) then
                        part11=smoothvalue1(i,k,0,0)
                        part12=smoothvalue2(i,k,0,0)
                    else
                        part11=smoothvalue1(i,k,itin,itout-1)
                        part12=smoothvalue2(i,k,itin,itout-1)
                    end if 
                    weight=(4.0*F(i,k)*(1.0-F(i,k))+1.0e-12)**0.25
                    do ii=-1,1
                        do kk=-1,1
                            sequence1(ii,kk)=(4.0*F(i+ii,k+kk)*(1.0-F(i+ii,k+kk))+1.0e-12)**0.25
                            sequence2(ii,kk)=sequence1(ii,kk)*smoothvalue1(i+ii,k+kk,itin,itout)
                            sequence3(ii,kk)=sequence2(ii,kk)*smoothvalue2(i+ii,k+kk,itin,itout)
                        end do 
                    end do 
                    call CNC_scheme(sequence1,part2down)
                    call CNC_scheme(sequence2,part21up)
                    call CNC_scheme(sequence3,part22up)
                    part21=part21up/part2down
                    part22=part22up/part2down
                    smoothvalue1(i,k,itin+1,itout)=weight*part11+(1.0-weight)*part21
                    smoothvalue2(i,k,itin+1,itout)=weight*part12+(1.0-weight)*part22
                    
                    Rcomponent=smoothvalue1(i,k,itin+1,itout)
                    Zcomponent=smoothvalue2(i,k,itin+1,itout)
                    mold=sqrt(Rcomponent**2+Zcomponent**2)
                    !re_normlization
                    smoothvalue1(i,k,itin+1,itout)=Rcomponent/(mold+1.0e-20)
                    smoothvalue2(i,k,itin+1,itout)=Zcomponent/(mold+1.0e-20)
                end do 
            end do 
            !extend to get the total smooth values
            do k=ZS,ZE
                i=RS-1
                smoothvalue1(i,k,itin+1,itout)=-smoothvalue1(i+1,k,itin+1,itout)
                smoothvalue2(i,k,itin+1,itout)=smoothvalue2(i+1,k,itin+1,itout)
                i=RE+1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i-1,k,itin+1,itout)
                smoothvalue2(i,k,itin+1,itout)=smoothvalue2(i-1,k,itin+1,itout)
            end do 

            do i=RS-1,RE+1
                k=ZS-1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i,k+1,itin+1,itout)
                smoothvalue2(i,k,itin+1,itout)=smoothvalue2(i,k+1,itin+1,itout)
                k=ZE+1
                smoothvalue1(i,k,itin+1,itout)=smoothvalue1(i,k-1,itin+1,itout)
                smoothvalue2(i,k,itin+1,itout)=smoothvalue2(i,k-1,itin+1,itout)
            end do 
            !update the values
            do i=0,NR
                do k=0,NZ
                    Nnormalcr(i,k)=smoothvalue1(i,k,itin+1,itout)
                    Nnormalcz(i,k)=smoothvalue2(i,k,itin+1,itout)
                end do 
            end do 
            
        end do 
    end do 

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

    !compute flux of every face

subroutine timeIntegratedFlux
    implicit none
    use par 
    !.....................................
    !R-direction face
    RS=1
    RE=NR
    
    ZS=1
    ZE=NZ-1
    do i=RS,RE
        do k=ZS,ZE
            !boundary face
            if(i==RS.or.i==NR) then
                Rvofflux(i,k)=0.
            else
                



subroutine isoadvector
    use par
    implicit none
    real(kind=dp)::dVf(4)
    
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1

    do i=RS,RE
        do k=ZS,ZE
            !initial dVf with upwind values
            if(U(i,k)>0.) then
                dVf(1)=F(i-1,k)*Phiw*dt
            else
                dVf(1)=F(i,k)*Phiw*dt
            end if 

            if(U(i+1,k)>0.) then
                dVf(2)=F(i,k)*Phie*dt
            else
                dVf(2)=F(i+1,k)*Phie*dt
            end if 

            if(W(i,k)>0.) then
                dVf(3)=F(i,k-1)*Phib*dt
            else
                dVf(3)=F(i,k)*Phib*dt
            end if 

            if(W(i,k+1)>0.) then
                dVf(4)=F(i,k)*Phit*dt
            else
                dVf(4)=F(i,k+1)*Phit*dt
            end if 
        end do 
    end do 

    do i=RS,RE
        do k=ZS,ZE
            if(cell_state(F(i,k)==0)) then
                
    return

end subroutine isoadvector

subroutine timeIntegratedFlux
    

subroutine limitFluxes

subroutine vof_update
    use par 
    implicit none
    RS=1
    RE=NR-1

    ZS=1
    ZE=NZ-1
    !initial vofflux with upwind value
    do i=RS,RE 
        do k=ZS,ZE
    !compute and update the vof_value
    do i=RS,RE
        do k=ZS,ZE
            Vp=RC(i)*dr*dz
            !interface cell
            if(cell_state(F(i,k))==0) then
                if(Phiwr(i,k)<0.) then !downwind face
                    call timeintegerflux(voffluxr(i,k))
                end if 
                if(Phiwr(i+1,k)>0.) then
                    call timeintegerflux(voffluxr(i+1,k))
                end if 
                if(Phiwz(i,k)<0.)
                    call timeintegerflux(voffluxz(i,k))
                end if 
                if(Phiwz(i,k)>0.)
                    call timeintegerflux(voffluxz(i,k+1))
                end if 
            end if 
            !boundary flux
            voffluxr(RS,k)=0.
            voffluxr(NR,k)=0.
            voffluxz(i,ZS)=0.
            voffluxz(i,NZ)=0.
            !get vof update
            dvf(1)=voffluxr(i,k)
            dvf(2)=voffluxr(i+1,k)
            dvf(3)=voffluxz(i,k)
            dvf(4)=voffluxz(i,k+1)
            F(i,k)=F(i,k)-(dvf(1)+dvf(2)+dvf(3)+dvf(4))/Vp

            if(F(i,k)<0..or.F(i,k)>1.) then
                call flux_correct
            end if 
        end do 
    end do 
