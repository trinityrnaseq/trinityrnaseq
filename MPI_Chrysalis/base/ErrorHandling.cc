#include "base/ErrorHandling.h"
//#include <malloc.h>


//#include <execinfo.h>

void print_trace(FILE *out, const char *file, int line)
{
  /*    const size_t max_depth = 100;
    size_t stack_depth;
    void *stack_addrs[max_depth];
    char **stack_strings;

    stack_depth = backtrace(stack_addrs, max_depth);
    stack_strings = backtrace_symbols(stack_addrs, stack_depth);

    fprintf(out, "Call stack from %s:%d:\n", file, line);

    for (size_t i = 1; i < stack_depth; i++) {
        fprintf(out, "    %s\n", stack_strings[i]);
    }
    free(stack_strings); // malloc()ed by backtrace_symbols
    fflush(out);*/
}


void SException::Throw() {
  print_trace(stdout, "", -1);

  throw *this;
}
  
