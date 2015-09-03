/*
 * Interview test for General Atomics-ASI
 * Nathan Honka
 * 6/8/2012
 *
 * The following code:
 *   - 
 */

 
 /****  Preprocessor ****/
 
// Includes
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
 
// Write file output parameters
#define WRITE_FILENAME 	"out.txt"
#define WRITE_FILE_TYPE	"w"
 
// Defined constants
#define DELTA_T					0.02	// s - discrete dT
#define PERIOD					4		// s - period of sine wave
#define TIME_LENGTH				20		// s - total length of input sine wave
#define T0						0		// s - signal start time
#define FILTER_STORAGE_LENGTH	3		// number of stored input & filtered values for filter
#define DIMENSIONS				3		// # of dimensions of the vector space in question

// Coefficients of G(z)
#define G_COEFF_A		0.0
#define G_COEFF_B		0.1218
#define G_COEFF_C		0.09338
#define G_COEFF_D		(-1.2382)
#define G_COEFF_E		0.4531
								
								
/**** Pre-allocated constants ****/
// Sine argument modifier equal to 2*pi/period
const double argMod = 2*M_PI/PERIOD;

// Navigation frame to body frame DCM
const double nav2body[DIMENSIONS][DIMENSIONS] =			
	{	{0.866, 0.25,	0.433	},
		{-0.5,	0.433,	0.75	},
		{0,		-0.866,	0.5		}	};

		
		
/**** Private function prototypes ****/

/*
 * Summary		Run the z-transform of G(z) = (a*z^2 + b*z + c) / (z^2 + d*z + e) on 
 * 				the input function u(t).
 * u			Pointer to array of length <arrayLength> containing input vector
 *				u(t).
 * y			Pointer to the array of length <arrayLength> containing the output
 *				vector y(t).
 * arrayLength	Integer equal to the length of the input vector u(t)
 * idx			Integer specifying the index of y to calculate via z-transform.
 */
static double zTransformG( double u[], double y[], int arrayLength );

/*
 * Summary		Transform a N-vector, multiplying by an NxN transformation matrix
 *				on the right.
 * vector		N-vector to be transformed.  Cannot be modified.
 * T			NxN transformation matrix.  Cannot be modified.
 * retVector	N-element array to store transformed vector.
 */
static void transformNVector( 	const double vector[DIMENSIONS],
								const double T[DIMENSIONS][DIMENSIONS],
								double retVector[DIMENSIONS] );


								
		
/****  Entry point  ****/

int main( void )
{
	// Allocate variables
	double u[FILTER_STORAGE_LENGTH];	// Storage for input sine wave with enough prior values for
										// z-transform.  Note that u[i] is U*z^(-i) (i.e., the ith previous value).
	double y[FILTER_STORAGE_LENGTH];	// Storage for filtered signal " " " ...
	double accNavVector[DIMENSIONS];	// Storage for acc vector, nav frame
	double accBodyVector[DIMENSIONS];	// Storage for acc vector, body frame
	double t;							// Current time (s)
	int i;								// Iterators
	FILE *pWriteFile;					// Text file write stream;

	
	/******* Initialization ********/
	// Initialize variables
	for( i=0; i<DIMENSIONS; i++ )
	{
		accNavVector[i] = 0;
		accBodyVector[i] = 0;
	}
	for( i=0; i<FILTER_STORAGE_LENGTH; i++ )
	{
		// Initialized to zero to fulfill causal filter requirement
		// ( u(t)=0 for t<0 )
		u[i] = 0;
		y[i] = 0;
	}
	t = 0;

	// Open write file stream & write file header
	pWriteFile = fopen (WRITE_FILENAME,WRITE_FILE_TYPE);
	fprintf( pWriteFile, "time\t\taccBodyX\t\t\taccBodyY\t\t\taccBodyZ\r\n" );
	
	
	/*** Create input waveform, filtered output, and transform to body coordinates. ***/
	/*** Print to file. ***/
	for( t=T0; t<=TIME_LENGTH; t+=DELTA_T )
	{	
		/****** Calculate input sine wave & filtered output ******/
		// Shift past storage arrays for z-transform
		for( i=FILTER_STORAGE_LENGTH-1; i>0; i-- )
		{
			u[i] = u[i-1];
			y[i] = y[i-1];
		}
		
		// Calculate input sine wave at current time t
		u[0] = sin(t*argMod);	// u = sin( 2*pi*t / period )

		// Calculate filtered output at current time t		
		y[0] = zTransformG( u, y, FILTER_STORAGE_LENGTH );
		
		/****** Calcualte body-frame acceleration ******/
		// Calculate nav-frame acceleration vector
		accNavVector[0] = y[0];
		accNavVector[1] = 0;
		accNavVector[2] = 0;
		
		// Transform acceleration vector to body frame
		transformNVector( accNavVector, nav2body, accBodyVector );
		
		/****** Print output to file ******/
		fprintf( pWriteFile, "%2.2f\t\t%2.10f\t\t%2.10f\t\t%2.10f\r\n",
			t, accBodyVector[0], accBodyVector[1], accBodyVector[2] );
	}
	
	// Close write stream
	fclose( pWriteFile );
	
	// Return
	return 0;
}




/****   Private functions   ****/


static double zTransformG( double u[], double y[], int arrayLength)
{
	double retVal = 0;

	// Notes on z-transform math:	
	// G(z) = (a*z^2 + b*z + c) / (z^2 + d*z + e)
	// 		= (a + b*z^-1 + c*z^-2) / (1 + d*z^-1 + e*z^-2 )
	// Difference equation: y[i] = sum( b_k * u[i-k] ) - sum( a_k * y[i-k] )
	// -> y[n] = a*u(n) + b*u(n-1) + c*u(n-2) - d*y(n-1) - e*y(n-2)
	// Note: causal filter: u[t]=0 for t<0.  This requirement is fulfilled by storage
	// arrays u & y that are initialized to 0.
	
	// Verify that array is greater than or equal to necessary time history for
	// z-transform.  In this case, since we require U*z^(-2), the array must have
	// at least the last 3 values.
	
	// If arrays are too small to perform z-transform in question, return 0.
	if( arrayLength < 3)
	{
		retVal = 0;
	}
	
	// Perform z-transform.
	else
	{
		retVal = G_COEFF_A*u[0] + G_COEFF_B*u[1] + G_COEFF_C*u[2]
			- G_COEFF_D*y[1] - G_COEFF_E*y[2];
	}
		
	return retVal;
}


static void transformNVector( 	const double vector[DIMENSIONS],
								const double T[DIMENSIONS][DIMENSIONS],
								double retVector[DIMENSIONS] )
{
	int		i,j;
	
	for( i=0; i<DIMENSIONS; i++ )
	{
		retVector[i] = 0;
		for( j=0; j<DIMENSIONS; j++ )
		{
			retVector[i] += vector[j]*T[i][j];
		}
	}

	return;
}