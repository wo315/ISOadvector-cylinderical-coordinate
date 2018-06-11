subroutine solver_possion
    use par
    implicit none
    !solve possion equation
    real(kind=dp),allocatable::A(:),B(:),X(:)
    integer,allocatable::JA(:),IA(:)
    integer::N,AI,NNR,NNZ,ik
    integer::iter,iprint,ijob,nhinv,nrest
    real(kind=dp)::tol
    !Íø¸ñ×ÜÊý
    NNR=NR-1
    NNZ=NZ-1
    N=NNR*NNZ
    !tolerance of residual pgmres
  
    tol=1.e-10_dp
    !amg
    ijob=0
    nrest=10
    
    nhinv=N
    iter=1000
    iprint=-6!only print the false information
    
    allocate(A(8*N),JA(8*N),IA(N+1),B(N),X(N))
    !allocate space

    
    AI=1
    IA(1)=1
    
    RS=1
    RE=NR-1
    ZS=1
    ZE=NZ-1
    !solve AX=B
    !AI·½³Ì×éÖÐÏµÊýµÄÎ»ÖÃ
    do k=ZS,ZE
            do i=RS,RE
                !ÏÈbottom£¬È»ºówest£¬Snorth,center,North,Top,±êºÅ´ÓÐ¡µ½´óÓë·½³ÌµÄÏµÊý¾ØÕóµÄÅÅÁÐÏà·ûºÏ
                ik=(k-ZS)*NNR+i-RS+1
                if(k>ZS) then! value at the bottom for it'is the smallest
                    A(AI)=-PAB(i,k)
                    JA(AI)=NAME(i,k-1)
                    AI=AI+1
                end if
                
                if(i>RS) then
                    A(AI)=-PAW(i,k)
                    JA(AI)=NAME(i-1,k)
                    AI=AI+1
                end if
            
                A(AI)=-PAP(i,k)
                JA(AI)=NAME(i,k)
                AI=AI+1
                
                if(i<RE) then
                    A(AI)=-PAE(i,k)
                    JA(AI)=NAME(i+1,k)
                    AI=AI+1
                end if
                
                if(k<ZE) then
                    A(AI)=-PAT(i,k)
                    JA(AI)=NAME(i,k+1)
                    AI=AI+1
                end if
                B(ik)=-PRHS(i,k)  !choose - for the equation
                
                IA(ik+1)=AI  !for the first IA(1)=1 ,IA(i+1)-IA(1) means the toal count number of AI
                X(ik)=0.  !innital values 
                end do
        end do
        
        ! write(*,*)' call amg solver to solve the poission equation'
        
        call DAGMG(N,A,JA,IA,B,X,ijob,iprint,nrest,iter,tol)
        
        RS=1
        RE=NR-1

        ZS=1
        ZE=NZ-1
        do k=ZS,ZE
                do i=RS,RE
                    ik=(k-ZS)*NNR+i-RS+1
                    Pm(i,k)=X(ik)
                end do
        end do
        
        deallocate(JA,IA,A,B,X)
        return
    end subroutine  solver_possion

subroutine RHS_possion
    use par
    implicit none
    
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    
    do i=RS,RE
            do k=ZS,ZE
                PRHS(i,k)=Um(i+1,k)*RP(i+1)*dz-Um(i,k)*RP(i)*dz+&
                Wm(i,k+1)*RC(i)*dr-Wm(i,k)*RC(i)*dr
                PRHS(i,k)=PRHS(i,k)/dt
            end do
    end do
    return 
    end subroutine RHS_possion


subroutine possion
    use par
    implicit none
    real(kind=dp)::denm
    real(kind=dp),pointer::Den(:,:)
    
    call RHS_possion
    
    !compute the cofficient of passion equation
    Den=>Denc
    RS=1
    RE=NR-1
    
    ZS=1
    ZE=NZ-1
    do i=RS,RE
            do k=ZS,ZE
                Denm=(Den(i,k)+Den(i,k)+Den(i,k+1)+Den(i,k+1))*0.25_dp !normal along with i
                PAW(i,k)=RP(i)*dz/Denm/dr
                
                Denm=(Den(i+1,k)+Den(i+1,k)+Den(i+1,k+1)+Den(i+1,k+1))*0.25_dp
                PAE(i,k)=RP(i+1)*dz/Denm/dr
                
                Denm=(Den(i,k)+Den(i,k)+Den(i+1,k)+Den(i+1,k))*0.25_dp
                PAB(i,k)=RC(i)*dr/Denm/dz
                
                Denm=(Den(i,k+1)+Den(i,k+1)+Den(i+1,k+1)+Den(i+1,k+1))*0.25_dp
                PAT(i,k)=RC(i)*dr/Denm/dz
                
                PAP(i,k)=-(PAW(i,k)+PAE(i,k)+PAB(i,k)+PAT(i,k))
            end do 
    end do
    call boundary_possion
    return
    end subroutine possion
