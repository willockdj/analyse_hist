/*******************************************************/
/*** Moves the place pointer in a read in int line *****/
/*** to just before the field requested            *****/
/*** Dave Willock April 1997                       *****/
/*******************************************************/
#include <stdio.h>

int next_space( int *p_ichar, int start, int num_of_chars );

int next_none_space( int *p_ichar, int start, int num_of_chars );

int find_field(int *p_ichar, int field_desired, int start, int num_of_chars )

{
int iloop;
int *p_this_char;
int place;

p_this_char= p_ichar+start;
place= start;

for (iloop=0; iloop < field_desired-1; iloop++)
  {
    place= next_none_space( p_this_char, place, num_of_chars); 

    place= next_space( p_this_char, place, num_of_chars);
  
  }

return place;
}
