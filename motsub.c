#include "slv_c_utils.h"
#include <stdio.h>

adams_c_Motsub    Motsub;

/*
* Note:  
* Use mixed case names for the Adams subroutine names when using the C 
* style interface.  For the default subroutine name capitalize the first
* letter and have the remaining letters lower case; Gfosub for example. 
* Doing this insures that Adams Solver correctly distinguishes a C style 
* subroutine from Fortran and calls with the appropriate interface.
*  
*/

void Motsub(const struct sAdamsMotion* motion, double time, int iord, int iflag, double* value)
{
	value[0] = 1;//曲柄驱动角速度的输入，大小为1rad/s
}
