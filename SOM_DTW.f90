program SOM_DTW_window
    implicit none
    character(len=50) :: cfg_name,input_name,bmu_out,weights_out,info_out
    integer*4 rr,cc,h1,h2,kk,c0,c1,d1,hd,w,win!i
    real*4 cost,last_min
    integer*4 i,j,k,l,ne,ig,jg,ir,ihit,jhit,tmp_idx,counter
    integer*4 n,m,m2,ny,nx,ndist,gridtype,neightype,shrinktype,learntype,nepochs1,nepochs2,nepochs,tmax
    real*4 radius,rmin,rmin2,rmax,amin,amax,amin2,rnd,stress,dhit,d,alpha,a2
    real*4, ALLOCATABLE, DIMENSION (:,:,:) :: som
    real*4, ALLOCATABLE, DIMENSION (:,:) :: dat,mat
    real*4, ALLOCATABLE, DIMENSION (:) :: a,b
    integer*4, ALLOCATABLE, DIMENSION (:) :: idx_n,mtchx,mtchy
    integer*4, dimension(37270) :: tmp_idx_vec = 0


    CALL get_command_argument(1, cfg_name)
    cfg_name = trim(cfg_name)
    print*, cfg_name
    open(1, FILE=cfg_name)
    read(1,*) input_name !input file name
    input_name = trim(input_name)
    read(1,*) bmu_out !bmu out file name
    bmu_out = trim(bmu_out)
    read(1,*) weights_out !weig out file name
    weights_out = trim(weights_out)
    read(1,*) info_out !training info out file name
    info_out = trim(info_out)
    read(1,*) n !number of objects
    read(1,*) m !length first series
    read(1,*) ny !number of SOM rows
    read(1,*) nx !number of SOM columns
    read(1,*) neightype !neighborhood type: gaussian (0) or bubble (1)
    read(1,*) shrinktype !neighborhood shrinking: exp (0) or lin (1)
    read(1,*) learntype !learning rate decay: exp (0) or lin (1)
    read(1,*) nepochs1 !training epochs
    read(1,*) nepochs2 !fine tuning epochs
    read(1,*) rmin2 !neighborhood radius during fine tuning phase
    read(1,*) rmin !neighborhood radius min value
    read(1,*) rmax !neighborhood radius max value
    read(1,*) amin !neighborhood alpha min value
    read(1,*) amax !neighborhood alpha max value
    read(1,*) win !window parameter Sakoe-Chiba
    close(1)

    w = max(win, abs(m-m))
    allocate(mat(m+1,m+1))

    !initialize the distance matrix
    do i = 1, (m+1)
        do j = 1, (m+1)
            mat(i,j) = 9.9E+09
        end do
    end do
    !set first value to 0
    mat(1,1) = 0

    !generate data index and read input file
    allocate(dat(n,m),idx_n(n))
    do i=1,n
        idx_n(i)=i
    end do

    open(2, file=input_name)
    do i=1,n
        read(2,*) dat(i,:)
    end do
    close(2)

    allocate(a(m),b(m))

    !set values of ancillary variables
    nepochs = nepochs1 + nepochs2
    tmax = nepochs1
    amin2 = 0.01
    allocate(mtchx(n), mtchy(n))

    !weights initialized with random data objects
    allocate(som(ny,nx,m))
    do i = 1,ny
        do j = 1,nx
            call RANDOM_NUMBER(rnd)
            tmp_idx = int(rnd*n)+1
            print*, tmp_idx
            som(i,j,:) = dat(tmp_idx,:)
        end do
    end do

    open(unit=3, file=info_out)
    do ne=1,nepochs
        stress = 0
        counter = 1
	do i = 1,n
       	    tmp_idx_vec(i)=0
        end do
	
        10 do while(tmp_idx_vec(n)==0)
            call RANDOM_NUMBER(rnd)
            tmp_idx = int(n * rnd) + 1
            do i = 1, counter
                if (tmp_idx_vec(i)==tmp_idx) then
                    go to 10
                end if
            end do
            tmp_idx_vec(counter) = tmp_idx
            counter = counter +1
        end do

        do i=1,n
            ihit = 0
            jhit = 0
            dhit = 9.9e30
            b = dat(tmp_idx_vec(i), :)

            do ig=1,ny
                do jg=1,nx
                    a(:) = som(ig,jg,:)

                    !DTW windowed distance
                    !set the values in the window range to 0
                    do k = 2, (m+1)
                        do l = max(2,(k-w)), min(m+1, (k+w))
                            mat(k,l) = 0.0
                        end do
                    end do

                    d = 0
                    !compute the distance between a and b (the series) in the window range
                    do k = 2, (m+1)
                        do l = max(2,(k-w)), min(m+1,(k+w))
                            cost = (a(k-1) - b(l-1))* (a(k-1) - b(l-1))
                            last_min = min(mat(k-1,l), mat(k,l-1), mat(k-1,l-1))
                            mat(k,l) = cost+last_min
                        end do
                    end do
                    d = mat(m+1,m+1)

                    if(d<dhit)then
                        dhit = d
                        ihit = ig
                        jhit = jg
                    end if
                end do
            end do
            stress = stress + dhit

            if(ne <= nepochs1) then
                if(shrinktype==1) then
                    radius = (rmax - rmin) * (tmax - ne) / tmax + rmin
                else
                    radius = rmax * (rmin / rmax) ** (ne / tmax)
                end if
                ir = int(radius)

                if(learntype==1) then
                    alpha = (amax - amin) * (tmax - ne) / tmax + amin
                else
                    alpha = amax * (amin / amax) ** (ne / tmax)
                end if
            else
                ir = int(rmin2)
                alpha = amin2
            end if

			!check if the row is even
            if((ihit / 2 == int(ihit / 2))) then
                h2 = 0
            else
                h2 = 1
            end if
            do rr = (ihit - ir),(ihit + ir)
				!check if the row exists
                if(rr >= 1.AND.rr <= ny)then
                    kk = ir - int(abs(ihit - rr) / 2)
					!check if row distance from BMU row is even
                    if(rr / 2 == INT(rr / 2)) then
                        h1 = 0
                    else
                        h1 = 1
                    end if
                    c0 = jhit - kk + h1 * (1 - h2)
                    c1 = jhit + kk - (1 - h1) * h2
                    do cc = c0,c1
						!check if the column exists
                        if(cc >= 1.AND.cc <= nx) then
                            if(neightype == 1)then
                                hd = abs(int(ihit) - int(rr))
                                if(abs(int(jhit) - int(cc)) > (abs(int(ihit) - int(rr)) / 2)) then
                                    hd = hd + abs(int(jhit) - int(cc)) - int(abs(int(ihit) - int(rr)) / 2)
                                end if
                                if(int(abs(int(ihit) - int(rr))) / 2 /= int(abs(int(ihit) - int(rr)) / 2)) then
                                    if((int(ihit) / 2) /= int(int(ihit) / 2)) then
                                        if(int(jhit) > (int(cc) + int(abs(int(ihit) - int(rr)) / 2)))then
                                            hd = hd - 1
                                        end if
                                    else
                                        if(int(cc) > (int(jhit) + int(abs(int(ihit) - int(rr)) / 2))) then
                                            hd = hd - 1
                                        end if
                                    end if
                                end if
                                d1 = hd
                                a2 = max(alpha * exp(-d1 * d1 / (2 * radius)), amin)!radius * radius
                            else
                                a2 = alpha
                            end if

                            do k = 1,m
                                som(rr, cc, k) = som(rr, cc, k) + a2 * (dat(tmp_idx_vec(i), k) - som(rr, cc, k))
                            end do
                        end if
                    end do !cc
                END IF
            end do !rr

        end do !next pattern

        !average stress = Quantization Error
        stress = stress/n
        write(3,"(9999(g0))") ne, ",", ir, ",", alpha, ",", stress

    end do !next epoch
    close(3)

    !final BMU
    do i=1,n
        ihit=0
        jhit=0
        dhit=9.9e30
        b = dat(i, :)
        do ig=1,ny
            do jg=1,nx
                d = 0
                a(:) = som(ig,jg,:)
                d = 0
                do k=1,m
                    d = d + (a(k) - b(k)) * (a(k) - b(k))
                end do
                d = sqrt(d)
                if(d<dhit)then
                    dhit = d
                    ihit = ig
                    jhit = jg
                end if
            end do
        end do
        mtchx(i) = jhit
        mtchy(i) = ihit
    end do!next pattern

    !writes the winning cell coord for each pattern
    open(unit=4, file=bmu_out)
    write(4,"(9999(g0))") n
    write(4,"(9999(g0))") ny, ",",nx
    do i=1,n
        write(4,"(9999(g0))") idx_n(i),",",mtchy(i),",",mtchx(i)
    end do
    close(4)

    open(unit=5, file=weights_out)
    write(5,"(9999(g0))") ny, ",",nx, ",",m
    do i=1,ny
        do j=1,nx
            do k = 1,m-1
                write(5,"(9999(g0,','))",advance='no') som(i,j,k)
            end do
            write(5,"(9999(g0))") som(i,j,m)
        end do
    end do
    close(5)

end program SOM_DTW_window




