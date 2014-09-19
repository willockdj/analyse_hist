#include <sys/times.h>
#include <unistd.h>
float etime_me(float *p_array)
{
    struct tms timestruct;
    long iret, sec, sec1;
    static long tick = 0;
    if(!tick)tick = sysconf(_SC_CLK_TCK);
    iret = times(&timestruct);
    sec = timestruct.tms_utime/tick;
    sec1 =  timestruct.tms_utime - sec*tick;
    *p_array = (float) sec + (float) sec1 / (float) tick;
    sec = timestruct.tms_stime/tick;
    sec1 =  timestruct.tms_stime - sec*tick;
    *(p_array+1) = (float) sec + (float) sec1 / (float) tick;
    return *p_array + *(p_array+1);
}

/*  Alternative forms 

float etime_(array)
     float array[];
{
  float etime();
  return etime(array);
}
float ETIME(array)
     float array[];
{
  float etime();
  return etime(array);
}
*/

