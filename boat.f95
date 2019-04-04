subroutine boatflight(i_size, j_size, range, ratio, prob, coast, cap_1_array, pop_1_array)
  implicit none
  integer, intent(in)  :: i_size, j_size, range
  real, intent(in) :: cap_1_array(:,:), coast(:,:), prob, ratio
  real, intent(inout) :: pop_1_array(:,:)
  real :: moving_pop, rand
  integer :: x, y, i, j, k, l

  ! write(*,*) filter
  ! Every grid location
  do i = 0, i_size-1
    do j = 0, j_size-1
      ! For every grid location
      ! Figure out if it is coast
      ! Find nearby grid points that are also coast, with 0 pop
      ! With probability prob, move ratio number of population
      ! To random found gridpoint

      ! Only consider coastal points
      if (coast(i+1,j+1) == 1) then
        ! Only if ratio > 0.95
        ! write(*,*) i,j
        ! write(*,*) (pop_1_array(i+1,j+1)/cap_1_array(i+1, j+1))
        if (((pop_1_array(i+1,j+1)/cap_1_array(i+1, j+1)) > 0.95) .AND.&
                          pop_1_array(i+1, j+1)>0.1) then
          ! write(*,*) "mass"
          moving_pop = ratio * pop_1_array(i+1, j+1)
          do k = 0, range-1
            do l = 0, range-1
              y = 1 + i + k - (range-1)/2
              x = 1 + j + l - (range-1)/2
              ! Check if the considered grid position is valid
              if ((y>0 .AND. y<=i_size) .AND. (x>0 .AND. x<=j_size)) then
                ! Check if the considered tile is an empty coast tile
                if ((coast(y,x) == 1) .AND. (pop_1_array(y,x) == 0) .AND. (cap_1_array(y,x)>0)) then
                  ! write(*,*) x,y,moving_pop, cap_1_array(y,x), pop_1_array(y,x)
                  ! write(*,*) cap_1_array(i+1,j+1), pop_1_array(i+1,j+1)
                  ! write(*,*) "eeyo"
                  ! With probability prob
                  call random_number(rand)
                  ! write(*,*) rand
                  if (rand <= prob) then
                    ! write(*,*) "Hello                               Move1!"
                    ! write(*,*) ratio
                    pop_1_array(i+1, j+1) = pop_1_array(i+1, j+1) - moving_pop
                    pop_1_array(y,x) = moving_pop
                  end if
                end if
              end if
            end do
          end do
        end if
      end if
      ! write(*,*) pop_1_array(i+1,j+1)
    end do
  end do
  ! write(*,*) "SUM",sum(trade_1_del)
  ! write(*,*) "SUM CAP", sum(cap_1_array)
  ! write(*,*) cap_1_array
end subroutine boatflight

subroutine boatflight2(i_size, j_size, range, ratio, prob, coast, cap_1_array, pop_1_array, pop_2_array)
  implicit none
  integer, intent(in)  :: i_size, j_size, range
  real, intent(in) :: cap_1_array(:,:), coast(:,:), prob, ratio
  real, intent(inout) :: pop_1_array(:,:), pop_2_array(:,:)
  real :: moving_pop, rand
  integer :: x, y, i, j, k, l

  ! write(*,*) filter
  ! Every grid location
  do i = 0, i_size-1
    do j = 0, j_size-1
      ! For every grid location
      ! Figure out if it is coast
      ! Find nearby grid points that are also coast, with 0 pop
      ! With probability prob, move ratio number of population
      ! To random found gridpoint

      ! Only consider coastal points
      if (coast(i+1,j+1) == 1) then
        ! Only if ratio > 0.95
        ! write(*,*) i,j
        ! write(*,*) (pop_1_array(i+1,j+1)/cap_1_array(i+1, j+1))
        if ((pop_1_array(i+1,j+1)/cap_1_array(i+1, j+1)) > 0.9) then
          ! write(*,*) "mass"
          moving_pop = ratio * pop_1_array(i+1, j+1)
          do k = 0, range-1
            do l = 0, range-1
              y = 1 + i + k - (range-1)/2
              x = 1 + j + l - (range-1)/2
              ! Check if the considered grid position is valid
              if ((y>0 .AND. y<=i_size) .AND. (x>0 .AND. x<=j_size)) then
                ! Check if the considered tile is an empty coast tile
                if ((coast(y,x) == 1) .AND. (pop_1_array(y,x) == 0) .AND. (cap_1_array(y,x)>0).AND. (pop_2_array(y,x)==0)) then
                  ! write(*,*) "eeyo"
                  ! With probability prob
                  call random_number(rand)
                  if (rand <= prob) then
                    ! write(*,*) "Move2!"
                    ! write(*,*) ratio
                    pop_1_array(i+1, j+1) = pop_1_array(i+1, j+1) - moving_pop
                    pop_1_array(y,x) = moving_pop
                  end if
                end if
              end if
            end do
          end do
        end if
      end if
      ! write(*,*) pop_1_array(i+1,j+1)
    end do
  end do
  ! write(*,*) "SUM",sum(trade_1_del)
  ! write(*,*) "SUM CAP", sum(cap_1_array)
  ! write(*,*) cap_1_array
end subroutine boatflight2
