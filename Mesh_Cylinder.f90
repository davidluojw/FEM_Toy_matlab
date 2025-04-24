	implicit real*8 (a-h, o-z)

    character*80 title
	allocatable :: x(:,:), node(:,:), ibctype(:,:), bcvalue(:,:)

    open(6, file='input.dat', status='unknown')
	open(7, file='input.plt', status='unknown')

    write(*,*) ' Give number of divisions for the quarter cylinder in r and theta directions:'
    read(*,*) n2, n1

    pressure = 1.d0
    R_i      = 1.d0
    R_o      = 2.d0
	dtheta = 90.d0/real(n1)
	dr     = (R_o-R_i)/real(n2)
    ! m  = 2*n1*n2  ! Triangle element
    m = n1*n2       ! Quadrilateral element
	n  = (n1+1)*(n2+1)

	allocate(x(2,n), node(4,m), ibctype(n,2), bcvalue(n,2))
    r     = R_i
    do j=1,n2+1
    do i=1,n1+1
      ii    = (j-1)*(n1+1) + i
      theta = dtheta*real(i-1)
      x(1,ii) = r*cosd(theta)
      x(2,ii) = r*sind(theta)
    enddo
    r     = r + dr
    enddo
    ! Triangle element
    ! do i=1,n1
    ! do j=1,n2
    !   kb = (j-1)*n1 + i
    !   i1 = (j-1)*(n1+1) + i
    !   i2 = i1 + 1
    !   i3 = i1 + (n1+1)
    !   i4 = i3 + 1
    !   k1 = 2*(kb-1) + 1
    !   node(1,k1) = i1
    !   node(2,k1) = i3
    !   node(3,k1) = i2
    !   k2 = k1 + 1
    !   node(1,k2) = i2
    !   node(2,k2) = i3
    !   node(3,k2) = i4
    ! enddo
    ! enddo

    ! Quadrilateral element
    do i=1,n1
    do j=1,n2
      kb = (j-1)*n1 + i
      i1 = (j-1)*(n1+1) + i
      i2 = i1 + 1
      i3 = i1 + (n1+1)
      i4 = i3 + 1
      node(1,kb) = i1
      node(2,kb) = i3
      node(3,kb) = i4
      node(4,kb) = i2
    enddo
    enddo

! Specify BC type and BC value:

    ibctype = 2
    bcvalue = 0.d0
    elem_len = dsqrt((x(1,1)-x(1,2))**2+(x(2,1)-x(2,2))**2)
    do i=1,n
      if(abs(dsqrt(x(1,i)**2+x(2,i)**2)-R_i) .le. 1.d-12) then
        costheta = x(1,i)/R_i
        sintheta = x(2,i)/R_i
        bcvalue(i,1) = pressure*elem_len*costheta
        bcvalue(i,2) = pressure*elem_len*sintheta
      endif
      if(i .eq. 1)    bcvalue(i,1) = bcvalue(i,1)/2.d0
      if(i .eq. n1+1) bcvalue(i,2) = bcvalue(i,2)/2.d0
      if(x(2,i) .eq. 0.d0)  then
        ibctype(i,2) = 1
      endif
      if(x(1,i) .eq. 0.d0) then
        ibctype(i,1) = 1
      endif
    enddo


! Output model:

	write(6, '(a20, i8, a10)') 	' A Cylinder Model with', m, ' Elements'
    write(6,*)  " Plane strain"
	write(6,*) m ,n,                    "              ! Numbers of elements and nodes"
    write(6,*)  "       1.d0       0.3d0               ! Young's modulus, Poisson's raio"
	write(6,*) '# Nodes:'
	do i=1,n
	  write(6,101) i, x(1,i), x(2,i)
	enddo
101 format(i10, 2d25.15)

    write(6,*) '# Element Connectivity:'
	do k=1,m
	  write(6,102) k, node(1,k), node(2,k), node(3,k), node(4,k)
	enddo
102 format(i10, 4i15)

    write(6,*) '# Load and Boundary Conditions:'
	do i=1,n
      write(6,103) i, ibctype(i,1), bcvalue(i,1), ibctype(i,2), bcvalue(i,2)
	enddo
103 format(i10, 2(5x,i10,5x,d25.15))

    write(6,*) '# End of File'

! Tecplot file:

	write(7,*) 'TITLE ="A Cylinder FEM Model"'
	write(7,*) 'VARIABLES = "X","Y" '
	write(7,*) 'ZONE T="2-D FE Domain", N=', n, ', E=', m, 'DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
	write(7,104) (x(1,i),x(2,i), i=1,n)
104 format(2E16.8)
	write(7,*)
	write(7,105) (node(1,k),node(2,k),node(3,k), k=1,m)
105 format(3I16)

	write(*,*) m, ' elements created!'

	end
