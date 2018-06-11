subroutine grid
    use par
    implicit none
    NR=NXD+1
    NZ=NZD+1
    allocate(RP(-2:NR+3),RC(-2:NR+3),ZP(-2:NZ+3),ZC(-2:NZ+3))

    Dr=ABS(XB-XA)/NXD
    Dz=ABS(YB-YA)/NZD

    RP(1)=XA
    RC(1)=RP(1)+dr*0.5_dp
    do i=2,NR+3
        RP(i)=RP(i-1)+dr
        RC(i)=RC(i-1)+dr
    end do
    !extend for velinterplote
    do i=0,-2,-1
        RP(i)=RP(i+1)-dr
        RC(i)=RC(i+1)-dr
    end do 

    ZP(1)=YA
    ZC(1)=ZP(1)+dz*0.5_dp
    do i=2,NZ+3
        ZP(i)=ZP(i-1)+dz
        ZC(i)=ZC(i-1)+dz
    end do

    do i=0,-2,-1
        ZP(i)=ZP(i+1)-dz
        ZC(i)=ZC(i+1)-dz
    end do 
    return
    end subroutine grid