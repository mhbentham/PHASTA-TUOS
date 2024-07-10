# CMake generated Testfile for 
# Source directory: /home/mli21/phasta-ncsu/phSolver/incompressible/test
# Build directory: /home/mli21/phasta-ncsu/phSolver/incompressible/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(clean_test_1 "./clean_tests.sh")
set_tests_properties(clean_test_1 PROPERTIES  WORKING_DIRECTORY "/home/mli21/phasta-ncsu/Test_Cases/")
add_test(submit_test "./submit_tests.sh")
set_tests_properties(submit_test PROPERTIES  WORKING_DIRECTORY "/home/mli21/phasta-ncsu/Test_Cases/")
add_test(Complete_Run_1 "python" "test.py")
set_tests_properties(Complete_Run_1 PROPERTIES  WORKING_DIRECTORY "/home/mli21/phasta-ncsu/Test_Cases/Shrey_Sample_Test/")
add_test(Complete_Run_2 "python" "test.py")
set_tests_properties(Complete_Run_2 PROPERTIES  WORKING_DIRECTORY "/home/mli21/phasta-ncsu/Test_Cases/test_bubble_tracking/")
add_test(Complete_Run_3 "python" "test.py")
set_tests_properties(Complete_Run_3 PROPERTIES  WORKING_DIRECTORY "/home/mli21/phasta-ncsu/Test_Cases/XYZTS_Sample_Test/")
