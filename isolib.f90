subroutine smooth_f
    use par
    implicit none
    integer::it,itm
    real(kind=dp)::Ffe,Ffw,Ffb,Fft
    real(kind=dp)::Afe,Afw,Afb,Aft
    real(kind=dp),pointer,dimension(:,:)::Ftem
    allocate(Ftem(1:NR-1,1:NZ-1))
    
    it=1
    itm=3
    
    do i=1,NR-1
            do k=1,NZ-1
                Ftem(i,k)=0_dp
            end do
    end do

    !......................................
    do while(it<itm)
        it=it+1 
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE
                    !the i,k respent the cell number
                    !get the smooth F
                    Ffw=0.5_dp*(FN(i,k)+FN(i-1,k))
                    Ffe=0.5_dp*(FN(i,k)+FN(i+1,k))
                    Ffb=0.5_dp*(FN(i,k)+FN(i,k-1))
                    Fft=0.5_dp*(FN(i,k)+FN(i,k+1))
                    
                    Afw=dz*RP(i)
                    Afe=dz*RP(i+1)
                    Afb=dr*RC(i)
                    Aft=dr*RC(i)

                    Ftem(i,k)=(Ffe*Afe+Ffw*Afw+Ffb*Afb+Fft*Aft)
                    Ftem(i,k)=Ftem(i,k)/(Afe+Afw+Afb+Aft)
                end do
        end do
        
        RS=1
        RE=NR-1
        ZS=1
        ZE=NZ-1 
        do i=RS,RE
                do k=ZS,ZE      
                    FN(i,k)=Ftem(i,k)
                end do
        end do
        
    end do
    
    !............................................................................ 
    !boundary condition use to compute the additional csf 
    call boundary_FN
    
    deallocate(Ftem)
    
    return
end subroutine smooth_f

subroutine order(A,n,B)
	!test ok //3_28
	implicit none
	integer,intent(in)::n
	integer,parameter::dp=selected_real_kind(p=15)
	real(kind=dp),dimension(n),intent(in)::A
	real(kind=dp),dimension(n),intent(out)::B
	integer::icy,jcy
	real(kind=dp)::temp
	!initial the values
	do icy=1,n
		B(icy)=A(icy)
	end do 
	! Sorting algorithm
	if(n==2) then
		if(B(1)>B(2)) then
			temp=B(1)
			B(1)=B(2)
			B(2)=temp
		end if 
	else
		do jcy=1,n-1
			!inner loop
			icy=1
			do while(icy<=n-jcy)
				if(B(icy)>B(icy+1)) then
					temp=B(icy)
					B(icy)=B(icy+1)
					B(icy+1)=temp
				end if 
			icy=icy+1
			end do 
			!out loop
		end do 
	end if 
	return 
end subroutine order

!solve the a*x=b matrix
!https://www.cfd-online.com/Wiki/LU_decomposition

subroutine LUsolver(a,b,n,x)
	!test ok2018_3_28
	implicit none
	integer,parameter::dp=selected_real_kind(p=15)
	integer,intent(in)::n
	real(kind=dp),dimension(n,n)::a
	real(kind=dp),dimension(n),intent(in)::b
	real(kind=dp),dimension(n),intent(out)::x
	real(kind=dp),dimension(n)::y
	integer::icy,jcy,kcy,scy
	real(kind=dp)::suma
	do jcy=1,n
		do icy=1,jcy
			if(icy==1) then
				a(icy,jcy)=a(icy,jcy)
			else
				suma=0.0
				do kcy=1,icy-1
					suma=suma+a(icy,kcy)*a(kcy,jcy)
				end do 
				a(icy,jcy)=a(icy,jcy)-suma
			end if 
		end do 
		if(jcy<n)	then
			do scy=1,n-jcy
				icy=jcy+scy
				if(jcy==1) then
					a(icy,jcy)=a(icy,jcy)/a(jcy,jcy)
				else
					suma=0.0
					do kcy=1,jcy-1
						suma=suma+a(icy,kcy)*a(kcy,jcy)
					end do 
					a(icy,jcy)=(a(icy,jcy)-suma)/a(jcy,jcy)	
				end if
			end do 
		end if 
	end do 

   !for now a(i,j) has been destory,stored L,U
	y(1)=b(1)
	do icy=2,n
		suma=0.0
		do jcy=1,icy-1
			suma=suma+a(icy,jcy)*y(jcy)
		end do
		y(icy)=b(icy)-suma 
	end do 

	x(n)=y(n)/a(n,n)
	do scy=1,n-1
		icy=n-scy
		suma=0.0
		do jcy=icy+1,n
			suma=suma+a(icy,jcy)*x(jcy)
		end do 
		x(icy)=(y(icy)-suma)/a(icy,icy)
	end do 

	return
end subroutine LUsolver 

 !cell_state==-1,vof=0
 !cell_state==0,interface cell
 !cell_state==1,full immersed

integer function cell_state(vofvalue,voftol)
	implicit none
	integer,parameter::dp=selected_real_kind(p=15)
	real(kind=dp)::voftol
	real(kind=dp)::vofvalue
	if(vofvalue<voftol) then
		cell_state=-1
	else if(vofvalue>1._dp-voftol) then
		cell_state=1
	else
		cell_state=0
	end if 
end function cell_state


subroutine velinterp(cell_label,isocenter,isoU0)
	!exp::inter to get the isocenter velocity
	use par 
	implicit none
	real(kind=dp),intent(in)::isocenter(2)
	integer,intent(in)::cell_label(2)
	real(kind=dp),intent(out)::isoU0(2)
	real(kind=dp)::w1,w2,weightsum,v1,v2,v3
	integer::a,b
	a=cell_label(1)
	b=cell_label(2)

	!interplote the isoU0(1) at the position(r,z) or position(x,z)
	w1=isocenter(1)-RP(a)
	w2=RP(a+1)-isocenter(1)
	weightsum=W1+w2
	if(weightsum==0._dp) then
		weightsum=weightsum+1.e-20
	end if 
	v1=(w2*U(a,b+1)+w1*U(a+1,b+1))/weightsum
	v2=(w2*U(a,b)+w1*U(a+1,b))/weightsum
	v3=(w2*U(a,b-1)+w1*U(a+1,b-1))/weightsum
	if(abs(isocenter(2)-Zc(b))<1.0e-20) then
		isoU0(1)=v2
	else if(isocenter(2)>Zc(b)) then
		w1=ZC(b+1)-isocenter(2)
		w2=isocenter(2)-Zc(b)
		isoU0(1)=(w2*v1+w1*v2)/(w1+w2)
	else if(isocenter(2)<Zc(b)) then 
	    w1=isocenter(2)-ZC(b-1)
		w2=Zc(b)-isocenter(2)
		isoU0(1)=(w2*v3+w1*v2)/(w1+w2)
	end if 

	!isoU0(2) interplote 
	w1=isocenter(2)-ZP(b)
	w2=ZP(b+1)-isocenter(2)
	weightsum=W1+w2
	if(weightsum==0._dp) then
		weightsum=weightsum+1.e-20
	end if 
	v1=(w2*W(a+1,b)+w1*W(a+1,b+1))/weightsum
	v2=(w2*W(a,b)+w1*W(a,b+1))/weightsum
	v3=(w2*W(a-1,b)+w1*W(a-1,b+1))/weightsum
	if(abs(isocenter(1)-RC(a))<1.0e-20) then
		isoU0(2)=v2
	else if(isocenter(1)>RC(a)) then
		w1=RC(a+1)-isocenter(1)
		w2=isocenter(1)-RC(a)
		isoU0(2)=(w2*v1+w1*v2)/(w1+w2)
	else if(isocenter(1)<RC(a)) then 
		w1=isocenter(1)-RC(a-1)
		w2=RC(a)-isocenter(1)
		isoU0(2)=(w2*v3+w1*v2)/(w1+w2)
	end if 

	return
end subroutine velinterp


subroutine nodesvof(label,nodesvalue)
	use par 
	implicit none
	integer,intent(in)::label(2)
	real(kind=dp),intent(out)::nodesvalue(4)

	nodesvalue(1)=Fnode(label(1),label(2))
	nodesvalue(2)=Fnode(label(1),label(2)+1)
	nodesvalue(3)=Fnode(label(1)+1,label(2)+1)
	nodesvalue(4)=Fnode(label(1)+1,label(2))

	return
end subroutine nodesvof

subroutine vertexloc(label,verloc)
	use par 
	implicit none
	integer,intent(in)::label(2)
	real(kind=dp),intent(out)::verloc(4,2)

	verloc(1,1)=RP(label(1))
	verloc(1,2)=ZP(label(2))

	verloc(2,1)=RP(label(1))
	verloc(2,2)=ZP(label(2)+1)

	verloc(3,1)=RP(label(1)+1)
	verloc(3,2)=ZP(label(2)+1)

	verloc(4,1)=RP(label(1)+1)
	verloc(4,2)=ZP(label(2))
	return
end subroutine vertexloc

subroutine cutpolygeo(inputisovalue,nodesvof,nodesloc,subcellapha)
	implicit none 
	integer,parameter::dp=selected_real_kind(p=15)
	integer::newmod,icy
	integer::polygeosize,itf,cut_number,number1,number2,casestate
	real(kind=dp)::polygeolist(5,2)
	real(kind=dp)::w1,w2,test,vp,deltar,deltaz
	real(kind=dp)::subdeltaz,subdeltar
	real(kind=dp),intent(in)::nodesloc(4,2)
	real(kind=dp),intent(in)::nodesvof(4)
	real(kind=dp),intent(in)::inputisovalue
	logical::face_state(4)
	real(kind=dp),intent(out)::subcellapha
    !get the cut_start
    cut_number=0
	do icy=1,4
		test=(inputisovalue-nodesvof(icy))*(inputisovalue-nodesvof(newmod(icy+1,4)))
		if(inputisovalue/=nodesvof(icy).and.test<=0.) then
			face_state(icy)=.true.
			cut_number=cut_number+1
			if(cut_number==1) then
				itf=icy
			end if 
		else
			face_state(icy)=.false.
		end if 
	end do 
	!start form it to get polygeo vertex
	polygeosize=0
	number1=newmod(itf,4)
	number2=newmod(itf+1,4)
	polygeosize=polygeosize+1
	w1=abs(inputisovalue-nodesvof(number1))
	w2=abs(inputisovalue-nodesvof(number2))
	polygeolist(polygeosize,1)=(nodesloc(number1,1)*w2+nodesloc(number2,1)*w1)/(w1+w2)
	polygeolist(polygeosize,2)=(nodesloc(number1,2)*w2+nodesloc(number2,2)*w1)/(w1+w2)
	if(nodesvof(number2)>inputisovalue) then
		casestate=1
	else
		casestate=0
	end if 
	polygeosize=polygeosize+1
	polygeolist(polygeosize,1)=nodesloc(number2,1)
	polygeolist(polygeosize,2)=nodesloc(number2,2)
	itf=itf+1
	!end of append the first two vertexs
	!get the end vertex as the second vertex
	cut_number=1
	do while(cut_number<2) 
		number1=newmod(itf,4)
		number2=newmod(itf+1,4)
		if(face_state(number1)) then
			polygeosize=polygeosize+1
			w1=abs(inputisovalue-nodesvof(number1))
			w2=abs(inputisovalue-nodesvof(number2))
			polygeolist(polygeosize,1)=(nodesloc(number1,1)*w2+nodesloc(number2,1)*w1)/(w1+w2)
			polygeolist(polygeosize,2)=(nodesloc(number1,2)*w2+nodesloc(number2,2)*w1)/(w1+w2)
			cut_number=cut_number+1
		else
			polygeosize=polygeosize+1
			polygeolist(polygeosize,1)=nodesloc(number2,1)
			polygeolist(polygeosize,2)=nodesloc(number2,2)
		end if 
		itf=itf+1
	end do !end of cut polygeo 
	!it compute the poly volume
	subcellapha=0.

	do icy=1,polygeosize
		number1=newmod(icy+1,polygeosize)!in this way it behave like a linke list
		subdeltaz=polygeolist(number1,2)-polygeolist(icy,2)!z2-z1
		subdeltar=(polygeolist(icy,1)**2+polygeolist(number1,1)**2+polygeolist(icy,1)*polygeolist(number1,1))
		subcellapha=subcellapha+1./6.*subdeltar*subdeltaz
	end do 
	
	if(casestate==1) then
		deltar=nodesloc(4,1)-nodesloc(1,1)  !r4-r1
		deltaz=nodesloc(2,2)-nodesloc(1,2)  !z2-z1
		vp=deltar*deltaz*(nodesloc(1,1)+nodesloc(4,1))*0.5_dp
		subcellapha=abs(subcellapha)/vp
	else if(casestate==0) then
		deltar=nodesloc(4,1)-nodesloc(1,1)  !r4-r1
		deltaz=nodesloc(2,2)-nodesloc(1,2)  !z2-z1
		vp=deltar*deltaz*(nodesloc(1,1)+nodesloc(4,1))*0.5_dp
		subcellapha=1.0-abs(subcellapha)/vp
	end if 
	return 
end subroutine cutpolygeo


subroutine isoinfo(inputisovalue,nodesvof,nodesloc,pos,centerpo,isonormal)
	implicit none
	!the defalut setting is that the heavy fuild is moving to the light fuild
	integer,parameter::dp=selected_real_kind(p=15)
	real(kind=dp),intent(in)::inputisovalue,nodesvof(4),nodesloc(4,2)
	real(kind=dp),intent(out)::pos(2,2),centerpo(2),isonormal(2)
	real(kind=dp)::test,w1,w2,Rcomponent,Zcomponent,mold,A,B
	logical::face_state(4)
	integer::itf,newmod,icy
	do icy=1,4
		test=(inputisovalue-nodesvof(icy))*(inputisovalue-nodesvof(newmod(icy+1,4)))
		if(inputisovalue/=nodesvof(icy).and.test<=0.) then
			face_state(icy)=.true.
		else
			face_state(icy)=.false.
		end if 
	end do 

	itf=0
	do icy=1,4
		if(face_state(icy)) then
			itf=itf+1
			w1=abs(inputisovalue-nodesvof(icy))
			w2=abs(inputisovalue-nodesvof(newmod(icy+1,4)))
			pos(itf,1)=(nodesloc(icy,1)*w2+nodesloc(newmod(icy+1,4),1)*w1)/(w1+w2)
			pos(itf,2)=(nodesloc(icy,2)*w2+nodesloc(newmod(icy+1,4),2)*w1)/(w1+w2)
		end if 
	end do 

	!center position
	centerpo(1)=(pos(1,1)+pos(2,1))*0.5_dp
	centerpo(2)=(pos(1,2)+pos(2,2))*0.5_dp
	!normal vector,make sure that the vector is direct from liquild to gas
	A=nodesvof(1)+nodesvof(2)-(nodesvof(3)+nodesvof(4))
	B=nodesvof(1)+nodesvof(4)-(nodesvof(2)+nodesvof(3))
	Rcomponent=sign((pos(2,2)-pos(1,2)),A)!z2-z1
	Zcomponent=sign((pos(2,1)-pos(1,1)),B)!r2-r1
	mold=sqrt(Rcomponent**2+Zcomponent**2)
	if(mold<1.0e-12) then
		mold=mold+1.0e-12
	end if 
	isonormal(1)=Rcomponent/mold
	isonormal(2)=Zcomponent/mold
	return
end subroutine isoinfo


