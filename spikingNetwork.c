/*******************************************************************************/
/* spikingNetwork.c */
/* Spiking network model proposed by */
/* Xiao-Jing Wang (2002) Neuron 36, 955-968 */
/* Shin Kira*/
/*******************************************************************************/
// Comment 3
// Then I added one...
// I am adding a comment here.
// Include external libraries:
#include <stdio.shin>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Define macros for the preprocessor:
#define MAXCHR1 101
#define synDelay 25
#define Eall 1600  // Total number of pyramidal neurons
#define Iall 400   // Total number of interneurons
#define Nall 2000 // Total number of pyramidal neurons

// Begin main function:
int main(int argc, char *argv[]) {



	int timeStep = atoi(argv[1]);
	int waitForInput = atoi(argv[2]);
	double  dt = atof(argv[3]);      // [ms]  Time step

	int i,j; // general index
	int t; // time step counter
	// For Shin's Computer: 
	// char OutFile0[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_spikes.txt"; // output file
	// char OutFile1[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_Ie.txt"; // output file
	// char OutFile2[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_Ii.txt"; // output file
	// char OutFile3[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_V.txt"; // output file
	// char OutFile4[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_S.txt"; // output file
	// char OutFile5[]="C:\\Users\\Shin\\Documents\\MATLAB\\DynamicalSystem\\SN_sumW.txt"; // output file

	// For Nick's Computer:
	char OutFile0[]="SN_spikes.txt"; // output file
	char OutFile1[]="SN_Ie.txt"; // output file
	char OutFile2[]="SN_Ii.txt"; // output file
	char OutFile3[]="SN_V.txt"; // output file
	char OutFile4[]="SN_S.txt"; // output file
	char OutFile5[]="SN_sumW.txt"; // output file
	char chrb[MAXCHR1]; // character buffer
	FILE *Fout0, *Fout1, *Fout2, *Fout3, *Fout4, *Fout5; // file pointer

	/* Main Settings*/
	double f  = 0.15;    // Fraction of selective pyramidal neurons
	double input1;
	double input2;
	double sumW[synDelay][4]={0};
	double V[Nall]={0};
	double IeA[Nall];
	double IrA[Nall];
	double IrN[Nall];
	double IrG[Nall];
	double Isyn[Nall];
	double SeA[Nall];  // S for external AMPA
	double WSrA[Nall]; // Weighted sum of S for recurrent AMPA
	double x[Eall];
	double SN[Eall];
	double SNtemp[Eall];
	double WSN1,WSN2,WSN3;
	double WSN[Nall];  // Weighted sum of S for recurrent NMDA
	double WSG[Nall];  // Weighted sum of S for recurrent GABA
	double bg_input[Nall];
	double spike_ind[Nall] = {0}; // Initial refractory period indicator
	int spikes[Nall] = {0}; // Initial spikes
	int spikeM[synDelay][Nall] = {0}; // Memory of spikes for a duration of synaptic delay.

	// Defining known parameters and setting initial conditions
	double Cmp   = 500;    // [pF]  Membrane capacitance of pyramidal neurons  
	double Cmi   = 200;    // [pF]  Membrane capacitance of interneurons
	double GLp   = 25;     // [nS]  Leak conductance of pyramidal neurons
	double GLi   = 20;     // [nS]  Leak conductance of interneurons

	// External and recurrent AMPA conductances in pyramidal cells
	double GeAp = 2.1;    // [nS]  Conductance of external AMPA current in pyramidal neurons
	double GrAp = 0.05;   // [nS]  Conductance of recurrent AMPA current in pyramidal neuron

	// External and recurrent AMPA conductances in interneurons
	double GeAi = 1.62;    // [nS]  Conductance of external AMPA current in interneurons
	double GrAi = 0.05;   // [nS]  Conductance of recurrent AMPA current in interneurons

	// NMDA and GABA conductances in pyramidal cells and interneurosns
	double GNp  = 0.165;   // [nS]  Conductance of recurrent NMDA current in pyramidal neuron
	double GGp  = 1.5;     // [nS]  Conductance of recurrent GABA current in pyramidal neuron ORIGNAL VALUE: 1.3, 1.5 works well...
	double GNi  = 0.13;   // [nS]  Conductance of recurrent NMDA current in interneurons
	double GGi  = 1.0;   // [nS]  Conductance of recurrent GABA current in interneurons
	double Mg = 1;         // [mM]  Concentration of magnesium
	double VL = -70;       // [mV]  Resting potential
	double Vre = -55;      // [mV]  Reset potential
	double Vth = -50;      // [mV]  Threshold
	double Ve = 0;         // [mV]  Reversal potential of excitatory synapse
	double Vi = -70;       // [mV]  Reversal potential of inhibitory synapse
	double TA = 2;         // [ms]  Time constant of AMPA receptors
	double TNr = 2;        // [ms]  Time constant of rising phase of NMDA receptors 
	double TNd = 100;      // [ms]  Time constant of decaying phase of NMDA receptors 
	double alpha = 0.5;    // [ms^-1]  Used to calculate SN below
	double TG  = 5;        // [ms]  Time constant of GABA receptors  

	// Defined by Nick:
	int counter;

	// Other settings:
	double wp = 1.7;            // Potentiated synaptic weight (originally used in the paper)
	double w  = 1;
	double wg = 1;				// Weights for inhibitory synapses = 1 as in the original paper.
								// It can be changed to see the effect.

	double randn(double mu, double sig); // Declare random number generator function, so that it can be moved to EOF


	
	// Definitions that rely on previously defined variables:
	int E1 = f * Eall;  // The number of selective population of pyramidal cells in S1
	int E2 = f * Eall;  // The number of selective population of pyramidal cells in S2
	double wn = 1 - f*(wp - 1)/(1 - f); // Depressed synaptic weight (originally used in the paper)


	//srand(50);					// Always use the same random numbers in simulations.
	
	srand(time(NULL));		// Set the argumet to "time(Null)" to generate different random numbers.

	printf("\nMaximum random variable: %d", RAND_MAX);
	if (((Fout0 = fopen(OutFile0,"wt")) && (Fout1 = fopen(OutFile1,"wt")) && (Fout2 = fopen(OutFile2,"wt")) && (Fout3 = fopen(OutFile3,"wt")) && (Fout4 = fopen(OutFile4,"wt")) && (Fout5 = fopen(OutFile5,"wt")) )==0) {
		strcpy(chrb,"The output file cannot be opened.  Stop the process.");
		printf("\n%s",chrb);
		return 0;
	}

	/** Setting initial values **/
	// All synapses have a latency of 0.5 msec as in the paper, 
	// S are gating variables at the synapses.
	// Parameters in a static state were saved previously.
	// Here, they are retrieved and used in the simulation.
	for (i=0;i<Nall;i++){
		Isyn[i] = -100 * (double) rand() / (double)RAND_MAX;
		SeA[i]  = 5    * (double) rand() / (double) RAND_MAX;
		V[i]    = Vre - 5 + 10 * (double) rand() / (double) RAND_MAX;
		WSrA[i] = 5    * (double) rand() / (double) RAND_MAX;
		if (i<Eall){
			x[i]    = 0;
			SN[i]   = 0; //0.3  * (double) rand() / (double) RAND_MAX;
		}
	}
	for (i=0;i<synDelay;i++){
		for (j=0;j<Nall;j++){
			spikeM[i][j] = ((double) rand() / (double)RAND_MAX) < 10/1000*dt ? 1:0;
		}
	}
	for (i=0;i<Nall;i++){
		WSG[i] = 10 * (double) rand() / (double) RAND_MAX;
	}

	/**********************************************************/
	/*************** Main Time-Stepping Loop ******************/
	/**********************************************************/

	for (t=0;t<timeStep-1;t++){
		/*Reset parameters*/
		counter=0;
		WSN1 = 0;
		WSN2 = 0;
		WSN3 = 0;

		// sumW	and spike M are stored in matrix so that the events in 0.5ms ago (synaptic delay) 
		// can be recalled to update synaptic parameters.
		/* Shifting sumW by one row. (sumW[24][:] corresponds to events 0.05ms ago) */
		for (i=1;i<synDelay;i++){
			for (j=0;j<=3;j++){
				sumW[i][j] = sumW[i-1][j];
			}
		}
		for (i=0;i<=3;i++){
			sumW[0][i]=0; // Reset the first row for new storage.
		}

		/* Shifting spikeM by one row. (spikeM[24][:] corresponds to events 0.05ms ago) */
		for (i=1;i<synDelay;i++){
			for (j=0;j<Nall;j++){
				spikeM[i][j] = spikeM[i-1][j];
			}
		}
		for (i=0;i<Nall;i++){
			spikeM[0][i]=0; // Reset the first row for new storage.
	}

		// Iterating membrane voltages and spikes
		for (i=0;i<Nall;i++){
			spike_ind[i] -= 1; // Neurons whose spike_ind is > 0 are in their refractory periods 
			if (spike_ind[i] < 0){
				spike_ind[i] = 0;  // If spike_ind goes below 0, set it to 0.
			}
			if (spike_ind[i] == 0 && i<Eall){
				V[i] += dt / Cmp * (-GLp*(V[i] -VL) - Isyn[i]); // Iterating the membrane voltage of pyramidal neurons
			}
			if (spike_ind[i] == 0){
				V[i] += dt / Cmi * (-GLi *(V[i]-VL) - Isyn[i]); // Iterating the membrane voltage of interneurons
			}
		}

		// Tracking membrane voltage of one representative neuron in each population
		fprintf(Fout3,"%5d%12.5f%12.5f%12.5f%12.5f\n",t,V[0],V[E1],V[E1+E2],V[Eall]); /* Tracking one neuron in each population. */

		for (i=0;i<Nall;i++){
			if ((V[i] > Vth) && (spike_ind[i] <= 0)){
				counter++;
				spikes[counter] = i;
				spikeM[0][i] = 1;

				// sumW is used for efficient calculations.
				// They are sum over the weights of neurons that spiked.
				if (i<E1){
					sumW[0][0] += wp;
					sumW[0][1] += wn;
					sumW[0][2] += w;
				}else if (i<(E1+E2)){
					sumW[0][0] += wn;
					sumW[0][1] += wp;
					sumW[0][2] += w;
				}else if (i<Eall){
					sumW[0][0] += wn;
					sumW[0][1] += wn;
					sumW[0][2] += w;	
				}else {
					sumW[0][3] += wg;
				}

				V[i] = Vre; // The voltage is reset to Vre = -55mV after the neuron spikes
				spike_ind[i] = (i <= Eall) ? 100: 50; // Implementing refractory periods (2ms = 100 steps).
				fprintf(Fout0,"%5d%5d\n",t,i); /* Print time and neuron# only when a neuron spikes */
			}

			// Additional external inputs (e.g. from the MT) to selective population 1 and 2.
			// The input value is rondomly drawn from a Gaussian distribution every 50ms (2500 steps). 
			if ((t % 2500)==0){
				input1 = (60 + randn(0,4))/1000; //[kHz]
				input2 = (20 + randn(0,4))/1000; //[kHz]
			}

			// No additional external imput before 100ms (5000 steps).
			if (t<waitForInput){
				bg_input[i] = (((double) rand() / (double) RAND_MAX) < (2.4*dt))? 1:0; // Background input 2400 Hz
			}else{
				if (i<E1){
					bg_input[i] = (((double) rand() / (double) RAND_MAX) < ((2.4 + input1)*dt))? 1:0; // 
				}else if (i<(E1+E2)){
					bg_input[i] = (((double) rand() / (double) RAND_MAX) < ((2.4 + input2)*dt))? 1:0; //
				}else{
					bg_input[i] = (((double) rand() / (double) RAND_MAX) < (2.4*dt))? 1:0; // Background input 2400 Hz
				}
			}

			SeA[i] += -SeA[i]*dt/TA + bg_input[i]; // Gating variables of external AMPA current

			/* Just printing the progress of iterations*/
			/*if (i==0){
				printf("\n t:%d",t);
			}*/
		}

		/** Update synaptic parameters with spikes that occured 0.5 ms ago. (sumW[24][:]) **/
		// WSA
		for (j=0;j<E1;j++){
			// WSrA[j] += -WSrA[j]*dt/TA + sumW[synDelay-1][0]; // Gating variables of recurrent AMPA current
			WSrA[j] += WSrA[j]*(-dt/TA + pow(-dt/TA,2)/2) + sumW[synDelay-1][0]; // Gating variables of recurrent AMPA current
		}
		for (j=E1;j<(E1+E2);j++){
			//WSrA[j] += -WSrA[j]*dt/TA + sumW[synDelay-1][1]; // Gating variables of recurrent AMPA current
			WSrA[j] += WSrA[j]*(-dt/TA + pow(-dt/TA,2)/2) + sumW[synDelay-1][1]; // Gating variables of recurrent AMPA current
		}
		for (j=(E1+E2);j<Nall;j++){
			//WSrA[j] += -WSrA[j]*dt/TA + sumW[synDelay-1][2]; // Gating variables of recurrent AMPA current
			WSrA[j] += WSrA[j]*(-dt/TA + pow(-dt/TA,2)/2) + sumW[synDelay-1][2]; // Gating variables of recurrent AMPA current
		}

		// WSG
		for (j=0;j<Nall;j++){
			//WSG[j]  += -WSG[j]*dt/TG + sumW[synDelay-1][3]; // Gating variables of recurrent GABA current
			WSG[j]  += WSG[j]*(-dt/TG + pow(-dt/TG,2)/2) + sumW[synDelay-1][3]; // Gating variables of recurrent GABA current
		}

		// WSN
		for (i=0;i<Eall;i++){
		//x[i] += -x[i]*dt/TNr  + spikeM[synDelay-1][i]; // Parameters used to calculate SN
		x[i] += x[i]*(-dt/TNr + pow(-dt/TNr,2)/2) + spikeM[synDelay-1][i]; // Parameters used to calculate SN
		}
		for (i=0;i<Eall;i++){
			//SN[i] += dt*(-SN[i]/TNd + alpha * x[i] * (1 - SN[i])); // Gating variables of recurrent NMDA current 
			SNtemp[i] = SN[i] + dt/2*(-SN[i]/TNd + alpha*x[i]*(1-SN[i]));
			SN[i] = SN[i] + dt*(-SNtemp[i]/TNd + alpha*x[i]*(1-SNtemp[i]));
		}

		/* WSN can be calculated in 3 patterns (see the synaptic weight matirx)*/
		for (i=0;i<Eall;i++){
			if (i<E1){
				WSN1 += wp*SN[i];
				WSN2 += wn*SN[i];
				WSN3 += w*SN[i];
			}else if (i<(E1+E2)){
				WSN1 += wn*SN[i];
				WSN2 += wp*SN[i];
				WSN3 += w*SN[i];	
			}else{
				WSN1 += wn*SN[i];
				WSN2 += wn*SN[i];
				WSN3 += w*SN[i];
			}
		}
		for (i=0;i<Nall;i++){
			if (i<E1)
				WSN[i] = WSN1;
			else if (i<(E1+E2))
				WSN[i] = WSN2;
			else
				WSN[i] = WSN3;
		}

		fprintf(Fout2,"%5d%12.5f%12.5f%12.5f\n",t,SN[0],SN[E1],SN[E1+E2]); // Check whether Sj are the same.
		fprintf(Fout4,"%5d%12.5f%12.5f%12.5f%12.5f\n",t,SeA[0],WSrA[0],WSN[0],WSG[0]); /* Tracking one neuron in each population. */
		fprintf(Fout5,"%5d%12.5f%12.5f%12.5f%12.5f\n",t,sumW[24][0],sumW[24][1],sumW[24][2],sumW[24][3]); /* Tracking one neuron in each population. */

		// Iterating synaptic currents
		for (i=0;i<Nall;i++){
			if (i<Eall){
				IeA[i] = GeAp * (V[i] - Ve) * SeA[i];
				IrA[i] = GrAp * (V[i] - Ve) * WSrA[i];
				IrN[i] = GNp * (V[i] - Ve) / (1 + Mg * exp(-0.062 * V[i]) / 3.57) * WSN[i];
				IrG[i] = GGp * (V[i] - Vi) * WSG[i];
			} else {
				IeA[i] = GeAi * (V[i] - Ve) * SeA[i]; // External AMPA current
				IrA[i] = GrAi * (V[i] - Ve) * WSrA[i]; // Recurrent AMPA current
				IrN[i] = GNi * (V[i] - Ve) / (1 + Mg * exp(-0.062 * V[i]) / 3.57) * WSN[i]; // Recurrent NMDA current
				IrG[i] = GGi * (V[i] - Vi) * WSG[i]; // Recurrent GABA current
			}
		}
		for (i=0;i<Nall;i++){
			Isyn[i] = IeA[i] + IrA[i] + IrN[i] + IrG[i]; // Net synaptic current
		}

		fprintf(Fout1,"%5d%15.5f%15.5f%15.5f%15.5f%15.5f\n",t,Isyn[0],IeA[0],IrA[0],IrN[0],IrG[0],V[0]); /* Tracking one neuron in each population. */

	}

	printf("\n Calculations completed.");
	fclose(Fout0);
	fclose(Fout1);
	fclose(Fout2);
	fclose(Fout3);
	return(0);
}

// Normal deviates generator (Box-Muller method)
// Spits out a random number whose mean is "mu" and standard deviation is "sig".
double randn(double mu, double sig)
{
	double v1, v2, rsq, fac;
	double storedval = 0;
	if (storedval == 0){
		do{
			v1 = 2 * (double) rand() / (double) RAND_MAX - 1;
			v2 = 2 * (double) rand() / (double) RAND_MAX - 1;
			rsq = v1 * v1 + v2 * v2;
		}while (rsq >= 1 || rsq == 0);
		fac = sqrt(-2 * log(rsq)/rsq);
		storedval = v1*fac;
		return (mu + sig*v2*fac);
	}else {
		fac = storedval;
		storedval = 0;
		return (mu + sig*fac);
	}
}
