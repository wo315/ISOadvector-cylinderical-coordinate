0,0.042,0.0,0.168             !# computational domain [XA,XB][YA,YB]
128,512                 !# number of nodes !even number
0.46,100				              !# terminal time,  number_write
3,1,1,1                                !# boundry condition 1=noslip 2=slip,3=axisymmetric,4=zero gradient,LEF,RIG,BOT,TOP
875.5, 0.118                           !# water:density  viscosity 1.0087e-3  
1.225, 1.77E-5                             !# air  :density  viscosity 1.81e-5 
1e-8                                   !# mass residual
0.0322                                  !# epiron
105,95                            	      !# thetaA,thetaR
3.36e-2,6.08e-3,7                                  !# h0,r0,deg
2e-5                             	  !#mindt
2                                        !NS_correction iteration max cycle number
1     !set_compute
