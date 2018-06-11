!interpolation method for advection term
!include upwind method and HR scheme
    subroutine HRscheme(Uf,Phi,PhiHR)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::Uf
    real(kind=dp),dimension(4),intent(in)::Phi
    real(kind=dp),intent(out)::PhiHR
    real(kind=dp)::PhiU,PhiD,PhiA
    real(kind=dp)::PhiDM,MM,KK
    
    if(Uf>=0._dp) then
        PhiU=Phi(1)
        PhiD=Phi(2)
        PhiA=Phi(3)
    else
        PhiU=Phi(4)
        PhiD=Phi(3)
        PhiA=Phi(2)
    end if
    
    if(abs(PhiA-PhiU)>1.0E-20) then
        PhiDM=(PhiD-PhiU)/(PhiA-PhiU)
    else
        PhiDM=0._dp
    end if
    
   !use the minmod scheme
   if (PhiDM>0._dp.and.PhiDM<=0.5_dp) then
       MM=1.5_dp
       KK=0_dp
   elseif(PhiDM>0.5_dp.and.PhiDM<=1._dp) then
       MM=0.5_dp
       KK=0.5_dp
   else
       MM=1._dp
       KK=0._dp
   end if
   PhiHR=MM*PhiD+KK*PhiA+(1._dp-MM-KK)*PhiU
   return
    end subroutine HRscheme
    
    subroutine upwind(Uf,Phi,PhiHR)
    implicit none
    integer,parameter::dp=selected_real_kind(p=15)
    real(kind=dp),intent(in)::Uf
    real(kind=dp),dimension(4),intent(in)::Phi
    real(kind=dp),intent(out)::PhiHR
    
    if(Uf>=0._dp) then
        PhiHR=Phi(2)!get the value of the upper cell
    else
        PhiHR=Phi(3)
    end if
    
    return
    end subroutine upwind
        
    