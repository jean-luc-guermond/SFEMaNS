MODULE COMM_LK_SFEM

        IMPLICIT NONE
        PUBLIC             :: LK_2_SFEM, SFEM_2_LK

    CONTAINS

        SUBROUTINE LK_2_SFEM(vec_sfem, vec_LK)
            USE def_type_field
            IMPLICIT NONE

            TYPE(mag_field_type), INTENT(IN)       :: vec_LK
            TYPE(mag_field_type), INTENT(OUT)      :: vec_sfem

            
            CALL vec_sfem%allocate_induction_fields
           
            vec_sfem%Hn = vec_LK%Hn
            vec_sfem%Bn = vec_LK%Bn
            vec_sfem%phin = vec_LK%phin
            vec_sfem%Hn1 = vec_LK%Hn1
            vec_sfem%Bn1 = vec_LK%Bn1
            vec_sfem%phin1 = vec_LK%phin1
            IF (.NOT. ALLOCATED(vec_sfem%time)) ALLOCATE(vec_sfem%time)
            vec_sfem%time = vec_LK%time
            vec_sfem%initialized = .TRUE.
        END SUBROUTINE LK_2_SFEM


        SUBROUTINE SFEM_2_LK(vec_LK, vec_sfem)
            USE def_type_field
            IMPLICIT NONE

            TYPE(mag_field_type), INTENT(IN)   :: vec_sfem
            TYPE(mag_field_type), INTENT(OUT)  :: vec_LK
             
            CALL vec_LK%allocate_induction_fields   

            vec_LK%Hn(:, :, :) = vec_sfem%Hn(:,:,:)
            vec_LK%Bn(:, :, :) = vec_sfem%Bn(:,:,:)
            vec_LK%phin(:, :, :) = vec_sfem%phin(:,:,:)
            vec_LK%Hn1(:, :, :) = vec_sfem%Hn1(:,:,:)
            vec_LK%Bn1(:, :, :) = vec_sfem%Bn1(:,:,:)
            vec_LK%phin1(:, :, :) = vec_sfem%phin1(:,:,:)
            ALLOCATE(vec_LK%time)
            vec_LK%time = vec_sfem%time
            vec_LK%initialized = .TRUE.


        END SUBROUTINE SFEM_2_LK


END MODULE COMM_LK_SFEM
