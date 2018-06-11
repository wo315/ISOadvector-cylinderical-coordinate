  subroutine next
    use par
    implicit none

    do i=0,NR+1
            do k=0,NZ+1
                Uo(i,k)=U(i,k)
                Wo(i,k)=W(i,k)
            end do
    end do

    do i=1,NR
        do k=1,NZ
            PO(i,k)=P(i,k)
        end do 
    end do 
    
    do i=-1,NR+1  !extend cells
            do k=-1,NZ+1
                FO(i,k)=F(i,k)
            end do
    end do

    return
    end subroutine next