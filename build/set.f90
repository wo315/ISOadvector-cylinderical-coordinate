subroutine set
    use par
    implicit none

    OPEN(1,FILE='./in.inc')

	READ(1,*) XA,XB,YA,YB
	READ(1,*) NXD,NZD
	READ(1,*) tend,number_write
	READ(1,*) BCXLEF,BCXRIG,BCYBOT,BCYTOP
	READ(1,*) Wdic,Wvis
	READ(1,*) Adic,Avis
	READ(1,*) mass
	READ(1,*) epiron
	READ(1,*) thetaA,thetaM
	READ(1,*) h0,r0,deg
	READ(1,*) mindt
	READ(1,*) itermaxzhu
	Read(1,*) set_compute
    

    close(1)

    return
    end subroutine set