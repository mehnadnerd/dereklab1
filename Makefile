test_mm: test_mm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG test_mm.c gen_matrix.c my_malloc.c -o test_mm

bmm: bmm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG -Werror -O3 -Ofast -ffast-math bmm.c gen_matrix.c my_malloc.c -o test_mm

asan: bmm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -fsanitize=address -g -DDEBUG -Werror -O0 bmm.c gen_matrix.c my_malloc.c -o test_mm

run_debug: bmm
	./test_mm 0 0 100

run_performance: bmm
	./test_mm 1 0 100

all:
	bmm

clean:
	rm *~; rm *.exe
