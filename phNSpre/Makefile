all: serial parallel

serial:
	cd src-serial; gmake 

parallel:
	cd src-parallel; gmake 

clean: clean-serial clean-parallel

clean-serial:
	cd src-serial; gmake clean

clean-parallel:
	cd src-parallel; gmake clean

rclean: rclean-serial rclean-parallel

rclean-serial:
	cd src-serial; gmake rclean

rclean-parallel:
	cd src-parallel; gmake rclean


