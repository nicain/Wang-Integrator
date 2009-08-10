######################################################### 
# PROGRAM NAME (no .f90) 
main=spikingNetwork
######################################################### 
CC=/usr/local/bin/gcc
objects1 = $(main).o 

$(main) : $(objects1) 
	$(CC) -lm -ftree-vectorize -o $(main) $(objects1) 

clean : 
	rm -f *.o $(main)
