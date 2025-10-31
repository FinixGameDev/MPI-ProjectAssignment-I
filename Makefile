all:
	mpicc main.c -lm -o fox
run:
	mpirun --use-hwthread-cpus -np 9 ./fox input1200
clean:
	rm -f fox output*
	rm -f output output*
	