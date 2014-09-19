/*********************************************************************************/
/**** Routine to make a triangular rep. matrix unit ******************************/
/**** Dave Willock April 1997                       ******************************/
/*********************************************************************************/

void set_unit_matrix(double *p_matrix, int length_side)
{

int index, side;

index= length_side;
 
while (index != 0)
  {
    *p_matrix= 1.0;
    for (side = index; side > 0; side--)
      {
        p_matrix++;
        *p_matrix= 0.0;
      }
    index--;
  }

return;
}
