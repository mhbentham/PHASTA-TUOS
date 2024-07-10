add_test(
        NAME clean_test_1
        COMMAND ./clean_tests.sh 
        WORKING_DIRECTORY ${PHASTA_SOURCE_DIR}/Test_Cases/
        )

add_test(
	NAME submit_test 
	COMMAND ./submit_tests.sh
	WORKING_DIRECTORY ${PHASTA_SOURCE_DIR}/Test_Cases/ 
	)

add_test(
	NAME Complete_Run_1
	COMMAND python test.py
	WORKING_DIRECTORY ${PHASTA_SOURCE_DIR}/Test_Cases/Shrey_Sample_Test/
	)	

add_test(
        NAME Complete_Run_2
        COMMAND python test.py
        WORKING_DIRECTORY ${PHASTA_SOURCE_DIR}/Test_Cases/test_bubble_tracking/
        )

add_test(
        NAME Complete_Run_3
        COMMAND python test.py
        WORKING_DIRECTORY ${PHASTA_SOURCE_DIR}/Test_Cases/XYZTS_Sample_Test/
        )


