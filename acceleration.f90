subroutine acceleration(nowtime,aga,agb)
    use par
    implicit none
    real(kind=dp),intent(in)::nowtime
    real(kind=dp),intent(out)::aga,agb
    aga=0._dp
    agb=-9.81_dp
  
    !agb=0._dp
    return
    end subroutine acceleration
    