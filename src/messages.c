#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "messages.h"

void message(const char *routine, const char *fmt, ...) {
  
  va_list ap;

  va_start(ap, fmt);  /* Point ap to first message. */
  
  printf("%s: ",routine);
  vprintf(fmt,ap);
  va_end(ap);
  printf("\n");
}


void errmess(const char *routine, const char *fmt, ...) {

  va_list ap;

  va_start(ap, fmt);
  printf("Error reported in: %s\n",routine);
  vprintf(fmt,ap);
  printf("\n Bombing out. \n");
  abort();
}

