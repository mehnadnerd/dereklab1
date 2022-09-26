#include<stdlib.h>
#include<inttypes.h>
#include "my_malloc.h"

int total_bytes_allocated = 0;

void *my_malloc(int size) {
    void *p;

    total_bytes_allocated += size;
    p = malloc(size + sizeof(uint64_t));  // assuming that malloc returns word aligned pointers
    ((uint64_t *) p)[0] = size;
    p += sizeof(uint64_t);
    return (p);
}

void *my_calloc(int num_elements, int element_size) {
    void *p;

    total_bytes_allocated += (num_elements * element_size);

    p = calloc(num_elements + 1, element_size);
    ((uint64_t *) p)[0] = num_elements * element_size;
    p += sizeof(uint64_t);
    return (p);
}

void my_free(void *p) {
    p -= sizeof(uint64_t);
    int size = ((uint64_t *) p)[0];
    total_bytes_allocated -= size;
    free(p);
}
