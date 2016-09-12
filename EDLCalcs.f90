subroutine calcvoltandfield(z_positions, z_plus, z_minus, eps, A_xy, valence, Voltage, Field, Bins, Nhalf)
    implicit none
    integer, intent(in) :: Bins, Nhalf, valence
    real(8), intent(in), dimension(0:Bins-1) :: z_positions
    real(8), intent(in), dimension(0:Nhalf-1) :: z_plus
    real(8), intent(in), dimension(0:Nhalf-1) :: z_minus
    real(8), intent(in) :: eps, A_xy
    real(8), intent(inout), dimension(0:Bins-1) :: Voltage
    real(8), intent(inout), dimension(0:Bins-1) :: Field
!f2py intent(in, out, inplace) :: Voltage,Field
    real(8) :: zbin
    integer :: i, j

    do i = 0, Bins - 1
	zbin = z_positions(i)
	do j=0, Nhalf-1
		if (z_plus(j)<=zbin) then
		Voltage(i) = Voltage(i) - (zbin - z_plus(j))*valence/(A_xy*eps)
		Field(i) = Field(i) - valence/(A_xy*eps)
		end if
	end do
    end do

    do i = 0, Bins - 1
	zbin = z_positions(i)
	do j=0, Nhalf-1
		if (z_minus(j)<=zbin) then
		Voltage(i) = Voltage(i) + (zbin - z_minus(j))*valence/(A_xy*eps)
		Field(i) = Field(i) + valence/(A_xy*eps)
		end if
	end do
    end do
end subroutine

subroutine Calc_EV(z_positions, Pos, exp_EV, exp_EV_HS, M, L_xy, epsWCA, sigWCA, Bins, NAtom)
    implicit none
    integer, intent(in) :: M, Bins, NAtom
    real(8), intent(in), dimension(0:Bins-1) :: z_positions
    real(8), intent(in), dimension(0:NAtom-1, 0:2) :: Pos  
    real(8), intent(in) :: epsWCA, sigWCA, L_xy
    real(8), intent(inout), dimension(1:Bins) :: exp_EV   !Note the different bin indexes
    real(8), intent(inout), dimension(1:Bins) :: exp_EV_HS
!f2py intent(in, out, inplace) :: exp_EV,exp_EV_HS
    real(8) :: zlo,zhi
    real(8) :: ranX,ranY,ranZ
    real(8), dimension(0:2) :: rij
    real(8) :: rij_scalar_sq,sigWCA_sq,rat6,sigHS,sigLJ_sq
    real(8) :: delU,mcount,delU_HS
    integer :: i, j
    call random_seed()
    sigWCA_sq = sigWCA**2
    sigHS = sigWCA*0.95
    sigLJ_sq = sigWCA_sq * 2**(-1./3)
    

    do i=1, Bins-1       
	zlo = z_positions(i-1)
	zhi = z_positions(i)

    	mcount = 0.
    	do while (mcount<M)
    		mcount = mcount + 1.
    		
		call random_number(ranX)
		call random_number(ranY)
		call random_number(ranZ)
		ranX = ranX*L_xy - 0.5*L_xy
		ranY = ranY*L_xy - 0.5*L_xy
		ranZ = ranZ*(zhi-zlo) + zlo  !This is correct!

    		delU = 0.
    		delU_HS = 0.
		do j=0,NAtom-1  !This loop pertains to only ONE insertion
			rij = Pos(j,:) - (/ ranX, ranY, ranZ /)		
			rij(:1) = rij(:1) - L_xy * dnint(rij(:1) / L_xy)
			rij_scalar_sq=sum(rij * rij)	
			if (rij_scalar_sq<=sigWCA_sq) then

				rat6 = (sigLJ_sq / rij_scalar_sq)**3 !This should be correct sigLJ_sq = sigWCA_sq * 2**(-1./3)
				
					! rat6 = (sigWCA_sq / rij_scalar_sq)**3 !This is the old way that incorrectly uses sigWCA, not sigLJ
				delU = delU + 4.*epsWCA*(rat6*(rat6-1)+0.25)
			end if

			rij_scalar_sq = rij_scalar_sq**0.5 !Did not want to name new variable
			if (rij_scalar_sq<=sigHS) then !There is at least one HS overlap
    				delU_HS = delU_HS + 1.
			end if
		end do
		exp_EV(i) = exp_EV(i) + exp(-delU)

		if (delU_HS.eq.0.) then  !There were no HS overlaps
			exp_EV_HS(i) = exp_EV_HS(i) + 1.	
		end if
	end do
	exp_EV(i) = exp_EV(i)/mcount
	exp_EV_HS(i) = exp_EV_HS(i)/mcount
    end do
end subroutine
