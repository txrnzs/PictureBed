/********************************************************************* 
UDF that simulates solidification by specifying a temperaturedependent viscosity property 
**********************************************************************/ 

#include "udf.h" 
DEFINE_PROPERTY(cell_viscosity1,c,t)
{ 
	real mu_lam; 
	real temp = C_T(c,t);
	if (temp > 488.) 
	{
		mu_lam = 5.5e-3; 
	}
	else if (temp > 486.) 
	{
		mu_lam = 143.2135 - 0.49725 * temp; 
	}
	else 
	{
		mu_lam = 1.; 
	}
	return mu_lam; 
}