program main
    use par
    implicit none
    call set
    call init
    call output(0)
    do  while(time<tend)
        call check_time
        call cal_Um
        call vel
        call check_ns
        call iso_vof
        call next
	    if(mod(nout,number_write)==0) THEN
	    call output(nout/number_write)
        call save_data
	    ENDIF
    end do
    stop
    end program main
