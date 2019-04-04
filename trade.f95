
subroutine calctrade(i_size, j_size, pop_1_array, cap_1_array, filter, trade_1_del, arguments)
  implicit none
  integer, intent(in)  :: i_size, j_size
  real, intent(in) :: pop_1_array(:,:), cap_1_array(:,:), filter(:,:)
  real, intent(inout) :: trade_1_del(:,:), arguments(:,:,:)
  real :: max, influence, local, foreign
  integer, dimension(2) :: args, filter_size
  integer :: x, y, i, j, k, l

  filter_size = shape(filter)
  ! write(*,*) filter
  ! Every grid location
  do i = 0, i_size-1
    do j = 0, j_size-1
      if (cap_1_array(i+1,j+1) >0) then
        ! Initialise at 0
        args = (/0, 0/)
        max = 0
        ! Find the largest influence in the filter
        ! If its population is larger than its own, transfer trade
        ! Every filter location
        do k = 0, filter_size(1)-1
          do l = 0, filter_size(1)-1
            ! Assign foreign coordinates
            y = 1 + i + k - (filter_size(1)-1)/2
            x = 1 + j + l - (filter_size(1)-1)/2

            ! Reset influence
            influence = 0
            ! Check if the considered location is withing the grid limits
            if ((y>0 .AND. y<=i_size) .AND. (x>0 .AND. x<=j_size)) then
              local = pop_1_array(1+i, 1+j)
              foreign = pop_1_array(y,x)
              ! Check if the foreign population is larger than the local one
              if (foreign>local) then
                influence = filter(k+1,l+1) * local * (1 - (local/foreign))
                ! influence = filter(y-i,x-j) * local * (1 - (local/foreign))
                ! Check if this is the largest influence
                if (influence > max) then
                  ! write(*,*) influence
                  max = influence
                  args(1) = y
                  args(2) = x
                end if
              end if
            end if
          end do
        end do
        ! Check if the trade transfer amount is larger than the cap
        ! if (cap_1_array(1+i,1+j) < 0) then
        !   write(*,*) cap_1_array(1+i,1+j)
        ! end if
        ! if (max>0) then
        !   write(*,*) max, cap_1_array(args(1), args(2))
        ! end if
        ! max = max * pop_1_array(1+i, 1+j)

        if (max>cap_1_array(1+i,1+j)) then
          ! write(*,*) max, cap_1_array(i+1,j+1)
          max = - cap_1_array(i+1,j+1)
        end if
        ! Once the greatest influencor is found, transfer trade
        trade_1_del(1+i,1+j) = trade_1_del(1+i,1+j) - max
        trade_1_del(args(1),args(2)) = trade_1_del(args(1),args(2)) + max
        arguments(1,1+i,1+j) = args(1)
        arguments(2,1+i,1+j) = args(2)
        ! if (max /= 0) then
          ! write(*,*) trade_1_del(1+i,1+j), trade_1_del(args(1),args(2))
          ! write(*,*) 1+i, 1+j, args, max
        ! end if
      end if
    end do
  end do
  ! write(*,*) "SUM",sum(trade_1_del)
  ! write(*,*) "SUM CAP", sum(cap_1_array)
  ! write(*,*) cap_1_array
end subroutine calctrade
