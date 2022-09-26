test_mm: test_mm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -g -DDEBUG -Werror -O3 -mno-mmx -mno-avx -mno-sse2 -fno-tree-vectorize test_mm.c gen_matrix.c my_malloc.c -o test_mm

asan: test_mm.c gen_matrix.c my_malloc.c gen_matrix.h my_malloc.h Makefile
	mpicc -fsanitize=address -g -DDEBUG -Werror -O0 test_mm.c gen_matrix.c my_malloc.c -o test_mm

run_debug: test_mm
	./test_mm 0 0 100

run_performance: test_mm
	./test_mm 1 0 100

all:
	test_mm

clean:
	rm *~; rm *.exe
