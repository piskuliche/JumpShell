Program JumpSolv
    implicit none
    integer :: i, j, k, ntimes, t, sk,cnt,to,n_tos
    integer :: acry_natms, nwater, nacryl
    real :: volume,maxdrsq, doh1, doh2
    integer, dimension(100) :: h2o_count, ad,h1, h2, hbond_count
    real, dimension (100) :: criteria,critsq
    character(len=3), dimension(100) ::  at_name
    character(len=2) :: ctmp2
    character(len=3) :: ctmp, ctmp3
    character(len=20) :: nfile
    integer, dimension(5000) :: hbondpartners,init_hbonds
    integer, dimension(5000000,5000) :: hbondlist
    real, dimension(500000,1000,3) :: eoh1, eoh2
    real, dimension(3) :: L, droh1, droh2
    real, dimension(1000,3) :: rO, r1, r2
    real, dimension(100,100,3) :: racryl
    real, dimension(5000) :: crp
    
    
    ! Set atomkey to 0.0  
    criteria=0.0
    rO=0.0; r1 = 0.0; r2 = 0.0
    racryl=0.0
    h2o_count=0
    critsq=0.0
    criteria=0.0
    h1 = 0; h2 = 0; ad = 0
    eoh1 = 0; eoh2 = 0
    hbondlist = 0; hbondpartners=0
    
    ! Reads the input file
    open(10,file='solv.in',status='old')
    read(10,*)
    read(10,*) nfile
    read(10,*) 
    read(10,*) volume, ntimes
    read(10,*)
    read(10,*) nwater, nacryl

    L(1)=volume ** (1.0/3.0)
    L(2)=L(1)
    L(3)=L(1)
    write(*,*) "Box length ", L(1)
    maxdrsq = 3*(L(1)/2)**2.0
    write(*,*) "Max Sq. Distance ", maxdrsq

    ! Reads the distances
    open(11,file="inp_vals", status='old')
    read(11,*) acry_natms
    write(*,*) "There are ", acry_natms, " per Acrylamide"
    cnt = 0
    do i=1,acry_natms
    read(11,*) at_name(i), criteria(i), ad(i), h1(i), h2(i)
    critsq(i)=criteria(i)**2.0
    if (criteria(i) .ne. 0.0) then
      cnt = cnt + 1
    endif
    enddo
    close(11)
    write(*,*) "There are ", cnt," atom criteria to compare"

    ! Opens the trajectory
    open(12,file=nfile, status='old')
    do t=1,ntimes
      ! Skip Header
      do sk=1,9
        read(12,*)
      enddo
      do i=1,nacryl
        do j=1,acry_natms
          read(12,*) ctmp, (racryl(j,i,k),k=1,3)
        enddo
      enddo
      do i=1,nwater
        ! Read Oxygen Atoms
        read(12,*) ctmp, (rO(i,k),k=1,3)
        ! Skip H Atoms (for now)
        read(12,*) ctmp, (r1(i,k),k=1,3)
        read(12,*) ctmp, (r2(i,k),k=1,3)
        r1(i,:) = r1(i,:) - L(:)*anint(( r1(i,:) - rO(i,:) )/L(:))
        r2(i,:) = r2(i,:) - L(:)*anint(( r2(i,:) - rO(i,:) )/L(:))
        droh1(:) = r1(i,:) - rO(i,:) - L(:)*anint(( r1(i,:) - rO(i,:) )/L(:))
        doh1   = sqrt(dot_product(droh1,droh1))
        eoh1(t,i,:)  = droh1(:)/doh1
        droh2(:) = r2(i,:) - rO(i,:) - L(:)*anint(( r2(i,:) - rO(i,:) )/L(:))
        doh2   = sqrt(dot_product(droh2,droh2))
        eoh2(t,i,:)  = droh2(:)/doh2
      enddo
      hbondpartners=0
      call CheckHbonds(nacryl,nwater,acry_natms, L, rO, r1,r2,racryl, h1, h2, hbondpartners)
      hbondlist(t,:)=hbondpartners(:)
    enddo! t loop
    close(12)

End Program   

Subroutine TCFLoop(ntos,sep,hbondlist,eoh1,eoh2)
  implicit none
  integer :: to, t, torigin
  integer :: i, j, k
  integer, dimension(5000) :: hbondpartners,init_hbonds,samehbond,laststep
  integer, dimension(5000000,5000) :: hbondlist
  real, dimension(5000) :: crp

  crp=0
  laststep=1
  do to=1,ntos
    torigin = to*sep
    init_hbonds(:)=hbondlist(torigin,:)
    do t=torigin,torigin+ncorr
      hbondpartners(:)=hbondlist(t,:)
      checkhbonds(init_hbonds,hbondpartners,laststep,samehbond)
      laststep(:)=samehbond(:) ! Sets the previous step
      ! Calculate indiv CRP
      do j=1,5000
        if (samehbond(j) .eq. 1) then
          crp(t) = crp(t) + 1
        endif
      enddo
      ! Need to calculate indiv P2(t) as well
      ! Need to sum stuff - should also consider whether I should
      ! Only calculate P2 if hbonds at time t (or at t0)
    enddo
  enddo

End Subroutine


Subroutine CheckPartners(init_partners,hbondpartners,laststep,samehbond)
  implicit none
  integer :: i, j, k, cnt
  integer, dimension(5000) :: init_partners,hbondpartners,samehbond,laststep
  samehbonds = 0 ! 0 if exchanged, 1 if same
  ! Note - laststep is exactly the same, just from the previous step
  do i=1,5000
    ! Check that to is not zero)
    if ( init_partners .eq. 0 .or. laststep(i) .eq. 0 ) then
      ! Doesn't count it - not hbonded at t=0 or if TCF is already 0
    else if ( hbondpartners .eq. 0 ) then
      ! Counts it as unbroken - transient break
      samehbonds(i)=1
    else
      ! Potential H-Bond Exchange - is it the same?
      if ( hbondpartners(i) .eq. init_partners(i)) then
        samehbonds(i)=1 ! Hbond Exchange
      endif
    endif
  enddo
   

End Subroutine


Subroutine CheckHbonds(nacryl, nwater, acry_natms, L, rO, r1, r2, racryl, h1, h2,critsq, hbondpartners)
  implicit none
  integer :: i, j, k, cnt, presenthbnd,acryindex
  integer :: nacryl, nwater, acry_natms
  integer :: oh1, oh2
  integer, dimension(100) :: h1, h2
  integer, dimension(5000) :: hbondpartners
  integer, dimension(1000) :: acryhbond
  real :: ang, dOO, dOx, dHX1, dHX2, dHO1, dHO2
  real :: rOXmax, rOOmax, rHOmax, rHXmax, angmax
  real, dimension(3) :: L
  real, dimension(100,100,3) :: racryl
  real, dimension (100) :: critsq
  real, dimension(1000,3) :: rO, r1, r2
  real, dimension(3) :: drOX, drHX1, drHX2,drHO1, drHO2, drOO
  real, dimension(3) ::eox, ehx1, ehx2, eho1, eho2, eoo
  ! This code calculates hbonding and stores it
  ! It stores this information in hbondpartners
  ! It calculates all the hydrogen bonds in the system (water-water,
  ! acryl->water, and water->acryl
  ! Note - acryhbond has nwater elements:
  !   These are 1 if hbonded to an acryl
  !   or 0 otherwise


  hbondpartners=0 
  rOXmax = 3.5; rHXmax = 2.45; angmax =30
  rOOmax = 3.1; rHOmax = 2.1
  hbondpartners=0
  acryhbond=0
  do k=1,nwater
    ! Checks Donations to Water j
    oh1=0;oh2=0
    do j=1,nwater
      if ( k .ne. j ) then
        if (oh1 .eq. 0 .or. oh2 .eq. 0) then
          drOO(:) = rO(k,:)-rO(j,:) - L(:)*anint((rO(k,:)-rO(j,:))/L(:))
          dOO = sqrt(dot_product(drOO,drOO)) 
          if  ( dOO < rOOmax ) then
            drHO1(:)=r1(k,:)-rO(j,:) - L(:)*anint((r1(k,:)-rO(j,:))/L(:)) 
            drHO2(:)=r2(k,:)-rO(j,:) - L(:)*anint((r2(k,:)-rO(j,:))/L(:))
            dHO1=sqrt(dot_product(drHO1,drHO1))
            dHO2=sqrt(dot_product(drHO2,drHO2))
            eoo(:) = drOO(:)/dOO
            if (dHO1 < rHOmax) then
              eho1(:)=drHO1(:)/dHO1
              ang = acosd(dot_product(eox,ehx1))
              if ( ang < angmax ) then
                oh1=j
              endif
            endif
            if (dHO2 < rHOmax) then
              eho2(:)=drHO2(:)/dHO2
              ang = acosd(dot_product(eoo,eho1))
              if (ang < angmax) then
                oh2=j
              endif
            endif ! dhO2 < rhomax
          endif ! oo < roomax
        endif ! oh1 oh2
      endif ! k ne j
    enddo ! j loop
    ! Checks Donations to Acryl j
    do i=1,acry_natms
      if ( critsq(i) .ne. 0.0 ) then
        do j=1, nacryl
          ! This index starts past the number of OHs available
          acryindex=j*acry_natms + i + nwater*2
          drOX(:)=rO(k,:)-racryl(i,j,:) - L(:)*anint((rO(k,:)-racryl(i,j,:))/L(:))
          dOX = sqrt(dot_product(drOX,drOX))
          if  ( dOX < rOXmax ) then
            ! Checks Water Donations to Acryl
            drHX1(:)=r1(k,:)-racryl(i,j,:) - L(:)*anint((r1(k,:)-racryl(i,j,:))/L(:))
            drHX2(:)=r2(k,:)-racryl(i,j,:) - L(:)*anint((r2(k,:)-racryl(i,j,:))/L(:))
            ! Water H, Acryl X (O or N) distance
            dHX1=sqrt(dot_product(drHX1,drHX1))
            dHX2=sqrt(dot_product(drHX2,drHX2))
            eox(:) = drOX(:)/dOX
            if ( dHX1 .lt. rHXmax ) then
              ehx1(:) = drHX1(:)/dHX1
              ang = acosd(dot_product(eox,ehx1))
              if ( ang < angmax ) then
                oh1 = acryindex
                acryhbond(k)=1 ! Water molecule is hbonded to acry
              endif
            endif
            if ( dHX2 .lt. rHXmax ) then
              ehx2(:) = drHX2(:)/dHX2
              ang = acosd(dot_product(eox,ehx2))
              if ( ang < angmax ) then
                oh2 = acryindex
                acryhbond(k)=1 ! Water molecule is hbonded to acry
              endif
            endif
          endif !  ( dOX < rOXmax )
          ! Checks Acryl Donations to Water
          if ( h1(i) .ne. 0) then
            drHO1(:)=racryl(h1(i),j,:) - rO(k,:) - L(:)*anint((racryl(h1(i),j,:) -rO(k,:))/L(:))
            ! Water H, Acryl X (O or N) distance
            dHO1=sqrt(dot_product(drHO1,drHO1))
            if ( dHO1 .lt. rHXmax ) then
              eho1(:) = drHO1(:)/dHO1
              ang = acosd(dot_product(-eox,eho1))
              if ( ang < angmax ) then
                hbondpartners(acryindex)=k
                acryhbond(k)=1 ! Water molecule is hbonded toacry
              endif
            endif
          endif
          if ( h2(i) .ne. 0) then
            drHO2(:)=racryl(h2(i),j,:) - rO(k,:) - L(:)*anint((racryl(h2(i),j,:)-rO(k,:))/L(:))
            dHO2=sqrt(dot_product(drHO2,drHO2))
            if ( dHO2 .lt. rHXmax ) then
              eho2(:) = drHO2(:)/dHO2
              ang = acosd(dot_product(-eox,eho2))
              if ( ang < angmax ) then
                hbondpartners(acryindex)=k
                acryhbond(k)=1 ! Water molecule is hbonded to acry
              endif
            endif
          endif
        enddo ! j nacryl
      endif ! crtisq
    enddo ! i acry_natms loop
    hbondpartners(2*k-1)=oh1
    hbondpartners(2*k)=oh2
  enddo ! k water loop

End Subroutine

