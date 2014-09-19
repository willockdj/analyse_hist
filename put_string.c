#include <stdio.h>
#include <limits.h>
#include "global_values.h"

/* routine to print out a line from a integer array */

void put_string( int *p_ichar, int length)
{

int *p_local;

p_local= p_ichar;

      while ( *p_local != '\0' && length != 0 )
      {
           putchar(*p_local);
           p_local++;
           length--;
      }

}

