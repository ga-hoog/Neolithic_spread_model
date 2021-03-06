
subroutine calcnj(i_size, j_size, pop_1_array, pop_2_array, filter, Nj_1, &
                  Nj_2, land)
  implicit none
  integer, intent(in)  :: i_size, j_size
  real, intent(in) :: pop_1_array(:,:), pop_2_array(:,:), filter(:,:), land(:,:)
  real, intent(inout) :: Nj_1(:,:), Nj_2(:,:)
  real :: old
  integer, dimension(2) :: filter_size
  integer :: x, y, i, j, k, l

  filter_size = shape(filter)


  ! Every grid location
  do i = 0, i_size-1
    do j = 0, j_size-1
      if (land(i+1, j+1) > 0) then
        ! Find Nj_y and Nj_x
        ! If its population is larger than its own, transfer trade
        ! Every filter location
        do k = 0, filter_size(1)-1
          do l = 0, filter_size(1)-1
            ! Assign foreign coordinates
            y = 1 + i + k - (filter_size(2)-1)/2
            x = 1 + j + l - (filter_size(2)-1)/2
            if ((y>0 .AND. y<=i_size) .AND. (x>0 .AND. x<=j_size)) then
              Nj_1(i+1,j+1) = Nj_1(i+1,j+1) + filter(k+1,l+1)*&
                                pop_1_array(i+1, j+1)*pop_2_array(y,x)
              Nj_2(i+1,j+1) = Nj_2(i+1,j+1) + filter(k+1,l+1)*&
                                pop_1_array(y, x)*pop_2_array(i+1,j+1)
            end if
          end do
        end do
      end if
    end do
  end do

end subroutine calcnj
