#include "slv_c_utils.h"
#include "stdio.h"
#include "userPortName.h"
#include "Windows.h"

adams_c_Consub    Consub;

/*
* Note:  
* Use mixed case names for the Adams subroutine names when using the C 
* style interface.  For the default subroutine name capitalize the first
* letter and have the remaining letters lower case; Gfosub for example. 
* Doing this insures that Adams Solver correctly distinguishes a C style 
* subroutine from Fortran and calls with the appropriate interface.
*  
*/

void Consub(REAL *par, int *npar)
{
	int i, j, nsize, errflg, nstates;
	int status;
	int iniflg = 1;
	int ipar[1];
	char errmsg[256];
	double states[3] = { 1,1,1 };
	FILE *outf;
	double time, dtime;
	time = 0;
	dtime = 0.06283;//每一步仿真0.06283s步长
	ipar[0] = 3;//摇杆的ID号
	outf = fopen("gwd.txt", "w");//摇杆的速度和加速度将输出到gwd.txt文件中
	fprintf(outf, "Wx(rad/s)\tWy(rad/s)\tWz(rad/s)\tAx(rad/s^2)\tAy(rad/s^2)\tAz(rad/s^2)\t\n");
	for (i = 0; i < 100; i++)//一共仿真100步
	{
		c_analys("KINEMATICS", "BY IGUO", time, time + dtime, iniflg, &status);//KINEMATICS运动学
		c_datout(&status);//完成一步仿真
		c_sysary("RVEL", ipar, 1, states, &nstates, &errflg);//提取摇杆的角速度信息
		for (j = 0; j < 3; j++)
			fprintf(outf, "%1f\t", states[j]);
		c_sysary("RACC", ipar, 1, states, &nstates, &errflg);//提取摇杆的角加速度信息
		for (j = 0; j < 3; j++)
			fprintf(outf, "%1f\t", states[j]);//将数据写入gwd.txt文件中
		fprintf(outf, "\n");
		time = time + dtime;

	}
	fclose(outf);

}
