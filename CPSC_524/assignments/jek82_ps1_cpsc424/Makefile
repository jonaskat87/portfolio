# Warning: Must run ''module load intel/2022b'' first.

CC = icc

#Various compiler options

CFLAG1 = -std=c99 
CFLAG2 = -g -O0 -fno-alias -std=c99 -fp-model precise
CFLAG3 = -g -O0 -fno-alias -std=c99
CFLAG4 = -g -O3 -no-vec -no-simd -fno-alias -std=c99
CFLAG5 = -g -O3 -xHost -fno-alias -std=c99

all: pi_1 pi_2 pi_3 pi_4 pi_5 div_1 div_2 div_3 div_4 div_5 triad

# Build code for Exercise 1
# Computing pi:
pi_1:	exercise_1.c timing.o
	$(CC) -o $@ $(CFLAG1) $^

pi_2:	exercise_1.c timing.o
	$(CC) -o $@ $(CFLAG2) $^

pi_3:	exercise_1.c timing.o
	$(CC) -o $@ $(CFLAG3) $^

pi_4:	exercise_1.c timing.o
	$(CC) -o $@ $(CFLAG4) $^

pi_5:	exercise_1.c timing.o
	$(CC) -o $@ $(CFLAG5) $^

# division
div_1:	test_div.c timing.o
	$(CC) -o $@ $(CFLAG1) $^

div_2:	test_div.c timing.o
	$(CC) -o $@ $(CFLAG2) $^

div_3:	test_div.c timing.o
	$(CC) -o $@ $(CFLAG3) $^

div_4:	test_div.c timing.o
	$(CC) -o $@ $(CFLAG4) $^

div_5:	test_div.c timing.o
	$(CC) -o $@ $(CFLAG5) $^

# Build code for Exercise 2 (uses CFLAG5, as "recommended for real codes")

triad:  exercise_2.c timing.o dummy.o
	$(CC) -o $@ $(CFLAG5) $^

# this removes executables and object files aside from timing.o and dummy.o
clean:
	rm -f pi* div* triad exercise*.o test_div.o
