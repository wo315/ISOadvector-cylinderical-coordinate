subroutine output(index_file)
  use par
  implicit none
  real(kind=dp)::X,Z,Uout,Wout,HH
  real(kind=dp)::masssum,massintx,massintz,centerx,centerz
  character(len=100)::filename
  integer::IO,index_file,MWW,Lfile

  PARAMETER(MWW=6)
  character*1 var_name(mww)*1,index_name*4
  character*8 title
  data var_name/"X","Z","U","W","P","vof"/
  write(index_name,'(I4.4)') index_file
  title="zvof"//index_name

  IO=100+cyc
  outtime=time
  masssum=0.
  massintx=0.
  massintz=0.
  centerx=0.
  centerz=0.
  write(*,*) index_file,'output tecplot files','*************************'

      open(10,position='Append',file='center.data')
      do i=1,NR-1
          do k=1,NZ-1
              massintx=massintx+RC(i)*dr*dz*(1.0-F(i,k))*RC(i)
              massintz=massintz+RC(i)*dr*dz*(1.0-F(i,k))*ZC(k)
              masssum=masssum+RC(i)*dr*dz*(1.0-F(i,k))
          end do
      end do
      centerx=massintx/masssum
      centerz=massintz/masssum
      write(10,100)masssum,centerx,centerz,time
100   format(4(f20.8)/)
      close(10)
      
      write(filename,"('file-'F8.4'.dat')")time

      open(unit=3,file=title//'.dat')
      write(3,*)"variables=",(',"',var_name(Lfile),'"',Lfile=1,MWW)
      write(3,110)NR-1,NZ-1,time
110   format('Zone i=',I3,',k=',I3,',F=POINT,solutiontime=',F8.4)
      do k=1,NZ-1
          do i=1,NR-1
              Uout=(U(i,k)+U(i+1,k))*0.5_dp
              Wout=(W(i,k)+W(i,k+1))*0.5_dp
              X=RC(i)
              Z=ZC(k)
              write(3,120)X,Z,Uout,Wout,P(i,k),F(i,k)
120           format(6(F20.8))
          end do
      end do
      close(3)
  return
  end subroutine output

  subroutine save_data
    use par 
    implicit none
    open(20,position='rewind',file='save.data')
    do k=1,NZ-1
            do i=1,NR
                write(20,*)U(i,k)
            end do    
        end do

    do i=1,NR-1
            do k=1,NZ 
                write(20,*)W(i,k)
            end do 
    end do 

    do i=1,NR-1
            do k=1,NZ-1
                write(20,*)P(i,k),F(i,k)
            end do    
    end do

    close(20)
    return 
    end subroutine save_data

    