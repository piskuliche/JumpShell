Program JumpSolv
    use, intrinsic :: iso_fortran_env
    implicit none
    integer :: i, j, k, ntimes, t, sk, cnt
    integer :: ncorr, to,ntos,totalconfigs, sep
    integer :: acry_natms, nwater, nacryl
    real :: volume,maxdrsq, doh1, doh2
    integer, dimension(100) :: h2o_count, ad,h1, h2, hbond_count
    real, dimension (100) :: criteria,critsq
    real, dimension (3) :: L
    character(len=3), dimension(100) ::  at_name
    character(len=2) :: ctmp2
    character(len=3) :: ctmp, ctmp3
    character(len=20) :: nfile

    
    write(*,*) "Begin Program: Jumps-Within-Shell"
    flush(OUTPUT_UNIT) 
    ! Set atomkey to 0.0  
    criteria=0.0
    h2o_count=0
    critsq=0.0
    criteria=0.0
    h1 = 0; h2 = 0; ad = 0
    ntos=100; sep=10; ncorr=2000; totalconfigs=100000
    
    ! Reads the input file
    open(10,file='solvshell.in',status='old')
    read(10,*)
    read(10,*) nfile
    read(10,*) 
    read(10,*) volume, ntimes
    read(10,*)
    read(10,*) nwater, nacryl
    read(10,*) 
    read(10,*) totalconfigs, sep, ncorr

    ntos=(totalconfigs-ncorr)/sep

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
    call TCFLoop(L, nacryl, acry_natms, nwater, critsq,ad,h1,h2, ntos,sep,ncorr)
    write(*,*) "End Program"
    close(12)

End Program   

Subroutine Read_Frame(L,nacryl, acry_natms, nwater,critsq,ad,h1,h2, hbondpartners, acryhbonds,e,siteloc)
  implicit none
  integer :: i, j, k
  integer :: sk, nacryl, acry_natms, nwater
  real    :: d1, d2
  character(len=3) :: ctmp
  
  integer, dimension(100) :: h1, h2, ad
  integer, dimension(5000) :: hbondpartners,acryhbonds,siteloc

  real, dimension (100) :: critsq

  real, dimension(1000,3) :: rO, r1, r2
  real, dimension(100,100,3) :: racryl

  real, dimension(1000, 3) :: e1,e2
  real, dimension(5000, 3) :: e
  real, dimension(3) :: L, dr1, dr2
  
  e1=0; e2=0; e=0
  siteloc=0
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
    dr1(:) = r1(i,:) - rO(i,:) - L(:)*anint(( r1(i,:) - rO(i,:) )/L(:))
    d1   = sqrt(dot_product(dr1,dr1))
    e1(i,:)  = dr1(:)/d1
    dr2(:) = r2(i,:) - rO(i,:) - L(:)*anint(( r2(i,:) - rO(i,:) )/L(:))
    d2   = sqrt(dot_product(dr2,dr2))
    e2(i,:)  = dr2(:)/d2
    ! Sets single vector with right dimensions to be sent out of code
    e(2*i-1,:) = e1(i,:)
    e(2*i,:) = e2(i,:)
  enddo
  hbondpartners=0;acryhbonds=0
  call CheckHbonds(nacryl,nwater,acry_natms, L, rO, r1,r2,racryl, h1, h2,critsq,ad, hbondpartners,acryhbonds,siteloc)

End Subroutine Read_Frame

Subroutine TCFLoop(L, nacryl, acry_natms, nwater, critsq,ad,h1,h2, ntos,sep,ncorr)
  implicit none
  integer :: to, t, torigin,ntos,sep,ncorr
  integer :: i, j, k, chkindex
  integer :: sk, nacryl, acry_natms, nwater
  integer :: count_crp
  integer :: count_acryc2, count_acryc2tos, count_acrycrp, count_acrycrptos
  real :: dp
  
  integer, dimension(5000) :: hbondpartners, acryhbonds, init_hbonds,samehbond,laststep,siteloc
  real, dimension(5000,3) :: e
  
  real, dimension(3) :: L
  integer, dimension(100) :: ad, h1, h2
  real, dimension (100) :: critsq

  integer, allocatable :: crptmp(:),acrycrptmp(:), bulkcrptmp(:)
  real, allocatable :: crp(:), acrycrp(:), c2(:), c2tmp(:)
  real, allocatable :: bulkcrp(:),bulkc2(:),bulkc2tmp(:)
  real, allocatable :: acryc2(:), acryc2tmp(:)
  integer, allocatable :: hbondlist(:,:), acryhbondlist(:,:), sitelist(:,:)
  real, allocatable :: elist(:,:,:)
  real, allocatable :: sitecrp(:,:), sitecrptmp(:,:)

  allocate( hbondlist(ncorr,5000) )
  allocate( acryhbondlist(ncorr,5000), sitelist(ncorr,5000) )
  allocate( crptmp(ncorr), acrycrptmp(ncorr), crp(ncorr), acrycrp(ncorr), c2tmp(ncorr), c2(ncorr) ) 
  allocate( acryc2(ncorr), acryc2tmp(ncorr) )
  allocate( elist(ncorr,5000,3) )
  allocate( sitecrp(ncorr,acry_natms), sitecrptmp(ncorr,acry_natms) )
  allocate( bulkcrp(ncorr), bulkcrptmp(ncorr), bulkc2(ncorr), bulkc2tmp(ncorr) )
  write(*,*) "Beginning TCF Loop"
  write(*,*) "Ncorr", ncorr
  write(*,*) "ntos", ntos
  write(*,*) "sep", sep
  ! Zero values that need to be zeroed
  hbondlist=0; acryhbondlist=0; elist=0.0;sitelist=0
  crp=0; acrycrp=0
  c2=0; acryc2=0
  bulkcrp=0; bulkc2=0
  sitecrp=0; siteloc=0; sitelist=0
  count_crp=0
  count_acryc2tos=0;  count_acrycrptos = 0
  
  ! Loop over time origins
  open(20,file="counts.dat")
  do to=1,ntos
    torigin = (to-1)*sep
    write(*,*) "Time Origin Reached:",to

    ! Zeros for time origin
    crptmp=0;acrycrptmp=0; c2tmp=0; acryc2tmp=0 ! These need to be zero every time
    bulkcrptmp=0; bulkc2tmp=0
    sitecrptmp=0
    laststep=1
    if (to .ne. 1) then ! Shifts by sep
      hbondlist = cshift(hbondlist,sep,dim=1)
      acryhbondlist = cshift(acryhbondlist,sep,dim=1)
      elist = cshift(elist,sep,dim=1)
      sitelist = cshift(sitelist,sep,dim=1)
    endif
    count_acryc2=0; count_acrycrp=0; count_crp=0
    do t=1,ncorr ! Loops over times
      if (to .eq. 1 .or. t .gt. ncorr-sep) then ! Either reads in the data
        call Read_Frame(L, nacryl, acry_natms, nwater, critsq,ad,h1,h2, hbondpartners, acryhbonds,e,siteloc)
        elist(t,:,:)=e(:,:)
        hbondlist(t,:)=hbondpartners(:)
        acryhbondlist(t,:)=acryhbonds(:)
        sitelist(t,:)=siteloc(:)
      endif
      call CheckPartners(hbondlist(1,:),hbondlist(t,:),laststep,samehbond)
      laststep(:)=samehbond(:) ! Sets the previous step

      ! Calculate indiv CRP + C2
      do j=1,5000
        if (j .le. nwater*2) then
          dp = dot_product(elist(t,j,:),elist(1,j,:))
          c2tmp(t) = c2tmp(t)+0.5*(3.0*dp**2.-1.0)
          acryc2tmp(t) = acryc2tmp(t) + 0.5*(3.0*dp**2.-1.0)*acryhbondlist(1,j)
          if ( t .eq. 1) count_acryc2=count_acryc2+acryhbondlist(1,j)
        endif
        if (samehbond(j) .eq. 1) then
          crptmp(t) = crptmp(t) + 1
          acrycrptmp(t) = acrycrptmp(t) + acryhbondlist(1,j)
          if ( t .eq. 1) then
            count_crp = count_crp + 1
            count_acrycrp = count_acrycrp + acryhbondlist(1,j)
          endif
          if (sitelist(t,j) .ne. 0) sitecrptmp(t,sitelist(t,j)) = sitecrptmp(t,sitelist(t,j)) + acryhbondlist(1,j)
        endif
      enddo
      ! Need to calculate indiv P2(t) as well
      ! Need to sum stuff - should also consider whether I should
      ! Only calculate P2 if hbonds at time t (or at t0)
    enddo
    ! Calc C2 Correlation 
    c2(:)=c2(:) + c2tmp(:)/(real(nwater)*2.0)
    acryc2(:)=acryc2(:) + acryc2tmp(:)/(real(count_acryc2))
    if ( acryc2tmp(1) .ne. 0.0 ) count_acryc2tos = count_acryc2tos + 1
    ! Calc CRP Correlation
    crp(:)=crp(:)+real(crptmp(:))/real(count_crp)
    acrycrp(:)=acrycrp(:)+real(acrycrptmp(:))/real(count_acrycrp)
    ! Calculate  bulk waters
    bulkc2tmp(:)  = c2tmp(:) - acryc2tmp(:)
    bulkcrptmp(:) = crptmp(:)-acrycrptmp(:)
    bulkc2(:) = bulkc2(:) + real(bulkc2tmp(:))/real(nwater*2 - count_acryc2)
    bulkcrp(:) = bulkcrp(:) + real(bulkcrptmp(:))/real(count_crp-count_acrycrp)
    if ( acrycrptmp(1) .ne. 0.0 )   count_acrycrptos = count_acrycrptos + 1
    write(20,*) t, nwater*2, count_acryc2, count_crp, count_acrycrp
    do k=1,acry_natms
      if (sitecrptmp(1,k) .ne. 0 ) then 
        sitecrp(:,k) = sitecrp(:,k) + real(sitecrptmp(:,k))/real(sitecrptmp(1,k))
      endif
    enddo
  enddo
  write(*,*) "Preparing to Dump Final Values"
  ! C2 Norm
  write(*,*) count_acryc2tos, acryc2(1), count_acrycrptos, acrycrp(1)
  c2(:)=c2(:)/real(ntos)
  acryc2(:)=acryc2(:)/real(count_acryc2tos)
  ! CRP Norm
  crp(:)=real(crp(:))/real(ntos)
  acrycrp(:)=real(acrycrp(:))/real(count_acrycrptos)
  bulkc2(:) =real(bulkc2(:))/real(ntos)
  bulkcrp(:)=real(bulkcrp(:))/real(ntos)
  do k=1,acry_natms
    if ( sitecrp(1,k) .ne. 0 ) sitecrp(:,k)=real(sitecrp(:,k))/real(sitecrp(1,k))
  enddo
  ! This writes to a file
  open(13,file="crp.dat")
  open(14,file="hbondcrp.dat")
  open(15,file="c2.dat")
  open(16,file="hbondc2.dat")
  open(17,file="sitecrp.dat")
  open(18,file="bulkcrp.dat")
  open(19,file="bulkc2.dat")
  do t=1,ncorr
    write(13,*) t, crp(t)
    write(14,*) t, acrycrp(t)
    write(15,*) t, c2(t) 
    write(16,*) t, acryc2(t)
    write(17,fmt='(100F12.5)') real(t), (sitecrp(t,k), k=1,acry_natms)
    write(18,*) t, bulkcrp(t)
    write(19,*) t, bulkc2(t)
  enddo
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  write(*,*) "End TCF Loop"

End Subroutine


Subroutine CheckPartners(init_partners,t_partners,laststep,samehbonds)
  implicit none
  integer :: i, j, k, cnt
  integer, dimension(5000) :: init_partners,t_partners,samehbonds,laststep
  ! This is a basic code to just check whether the h-bond is the same or
  ! different over the loop.

  samehbonds = 0 ! 0 if exchanged, 1 if same
  cnt = 0
  ! Note - laststep is exactly the same, just from the previous step
  do i=1,5000
    ! Check that to is not zero)
    if ( init_partners(i) .eq. 0 .or. laststep(i) .eq. 0 ) then
      ! Doesn't count it - not hbonded at t=0 or if TCF is already 0
    else if ( t_partners(i) .eq. 0 ) then
      ! Counts it as unbroken - transient break
      samehbonds(i)=1
    else
      ! Potential H-Bond Exchange - is it the same?
      if ( t_partners(i) .eq. init_partners(i)) then
        !write(*,*) 'Same Partner'
        samehbonds(i)=1 ! Hbond Exchange
      else
        cnt = cnt + 1
      endif
    endif
  enddo
  !write(*,*) cnt, "H-Bond Exchanges"
   

End Subroutine


Subroutine CheckHbonds(nacryl, nwater, acry_natms, L, rO, r1, r2, racryl, h1, h2,critsq,ad, hbondpartners,acryhbond,siteloc)
  implicit none
  integer :: i, j, k, cnt, presenthbnd,acryindex
  integer :: nacryl, nwater, acry_natms
  integer :: oh1, oh2, oh1index, oh2index
  integer :: acrydonorflag
  integer, dimension(100) :: h1, h2,ad
  integer, dimension(5000) :: hbondpartners
  integer, dimension(5000) :: acryhbond, siteloc
  real :: ang, dOO, dOx, dHX1, dHX2, dHO1, dHO2
  real :: rOXmax, rOOmax, rHOmax, rHXmax, angmax,anghohmax
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
  ! Note - acryhbond has the same dimensionality as hbondpartners:
  !   These are 1 if hbonded to an acryl
  !   or 0 otherwise

  hbondpartners=0
  acrydonorflag=0 ! 0 turns of considering acryl donations

  ! H-bond Criteria 
  rOXmax = 3.5; rHXmax = 2.45; angmax =30
  rOOmax = 3.1; rHOmax = 2.0; anghohmax = 20
  ! END H-bond Criteria


  hbondpartners=0
  acryhbond=0
  do k=1,nwater
    oh1index=2*k-1
    oh2index=2*k
    oh1=0;oh2=0
    ! Water-Water H-Bonds
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
              ang = acosd(dot_product(eoo,eho1))
              if ( ang < anghohmax ) then
                oh1=j
              endif
            endif
            if (dHO2 < rHOmax) then
              eho2(:)=drHO2(:)/dHO2
              ang = acosd(dot_product(eoo,eho2))
              if (ang < anghohmax) then
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
          acryindex=(j-1)*acry_natms + i + nwater*2
          drOX(:)=rO(k,:)-racryl(i,j,:) - L(:)*anint((rO(k,:)-racryl(i,j,:))/L(:))
          dOX = sqrt(dot_product(drOX,drOX))
          if  ( dOX < rOXmax .and. ad(i) .eq. 1) then
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
                acryhbond(oh1index)=1 ! Water molecule is hbonded to acry
                acryhbond(oh2index)=1 ! Water molecule is Hbonded to acry
                siteloc(oh1index)=i ! Stores which acryl atom nearby
                siteloc(oh2index)=i ! Stores which acryl atom nearby
              endif
            endif
            if ( dHX2 .lt. rHXmax ) then
              ehx2(:) = drHX2(:)/dHX2
              ang = acosd(dot_product(eox,ehx2))
              if ( ang < angmax ) then
                oh2 = acryindex
                acryhbond(oh1index)=1 ! Water molecule is hbonded to acry
                acryhbond(oh2index)=1 ! Water " ...  "
                siteloc(oh1index)=i ! Stores which acryl atom nearby
                siteloc(oh2index)=i ! Stores which acryl atom nearby
              endif
            endif
          endif !  ( dOX < rOXmax )
          ! Checks Acryl Donations to Water
          if (acrydonorflag .eq. 1) then
            if ( h1(i) .ne. 0) then
              drHO1(:)=racryl(h1(i),j,:) - rO(k,:) - L(:)*anint((racryl(h1(i),j,:) -rO(k,:))/L(:))
              ! Water H, Acryl X (O or N) distance
              dHO1=sqrt(dot_product(drHO1,drHO1))
              if ( dHO1 .lt. rHXmax ) then
                eho1(:) = drHO1(:)/dHO1
                ang = acosd(dot_product(-eox,eho1))
                if ( ang < angmax ) then
                  hbondpartners(acryindex)=k
                  acryhbond(acryindex)=1 ! Water molecule is hbonded toacry
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
                  acryhbond(acryindex)=1 ! Water molecule is hbonded to acry
                endif
              endif
            endif
          endif ! End Acryl donations to water
        enddo ! j nacryl
      endif ! crtisq
    enddo ! i acry_natms loop
    hbondpartners(2*k-1)=oh1
    hbondpartners(2*k)=oh2
  enddo ! k water loop

End Subroutine

