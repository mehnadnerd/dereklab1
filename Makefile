test_mm: test_mm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG test_mm.c gen_matrix.c my_malloc.c -o test_mm

bmm: bmm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG -Werror -O3 -Ofast -ffast-math bmm.c gen_matrix.c my_malloc.c -o test_mm

asan: bmm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -fsanitize=undefined -g -DDEBUG -Werror -O0 bmm.c gen_matrix.c my_malloc.c -o test_mm

cmm: cmm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG -Wall -Werror -O3 -Ofast -ffast-math -march=native cmm.c gen_matrix.c my_malloc.c -o test_mm

run_debug: cmm
	./test_mm 0 0 4

run_performance: cmm
	./test_mm 1 0 100

all:
	cmm

clean:
	rm *~; rm *.exe

dcilk:
	/scratch/04347/marui/cilk/build/bin/clang -fopencilk -fsanitize=cilk -Og -g -O3 cilk_mm.c gen_matrix.c my_malloc.c -o test_mm

cilk:
	/scratch/04347/marui/cilk/build/bin/clang -fopencilk -O3 cilk_mm.c gen_matrix.c my_malloc.c -o test_mm
