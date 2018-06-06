//Compile linking with fftw3:
// gcc <NAMEOFFILE>.c -lfftw3 -lm
//Execute:
// ./<NAMEOFFILE>.out filename

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

//*******************************************************************************************************************************************************


int main(int argc, char *argv[]){
  int i,j,k,l,T,n,cnt,nodo,numnodos,numsur,nbin = 128,dim,dim2,index,delay,**delaynodes,four=1,size;
  float *s,**phase_node,***phase_surr,*correltime,**correlnodes,*correlsurr,*correlsurrtime,hist[nbin],*correllow,*correlhigh,aux,*autocorr,pi,maxx,minn,dx,prob,threshold,threshold2;
  FILE *fp;
  char *line,*variables,*filename;
  size_t line_size=1024;
  fftw_complex *data1,*data2,*hlbrt1;
  fftw_plan pf,pb;

  //***********

  if (argc==2){
    size = strlen(argv[1]);
    filename = (char *) malloc(size * sizeof(char));
    strcpy(filename, argv[1]);
    printf("%s\n",filename);
  }
  else{
    printf("Usage: ./<executable_filename>.out <input_filename>");
    strcpy(filename, "out.csv");
  }
  
  //**********

  pi = acos(-1.);
  
  line = malloc(line_size * sizeof(char));
  variables = malloc(line_size * sizeof(char));

  //*********************** READING DATA PARAMETERS ********************************************
  
  fp=fopen(filename,"r");
  if (getline(&line, &line_size, fp) != -1){ //# channels = 67, surrogates = 100, lenght = 500
    if(line[0] == '#'){
      variables = strtok(line," \t");
      variables = strtok(NULL," \t");
      variables = strtok(NULL," \t");
      variables = strtok(NULL,",");
      numnodos = atoi(variables);
      printf("numnodos = %i\n", numnodos);
      
      variables = strtok(NULL," \t");
      variables = strtok(NULL," \t");
      variables = strtok(NULL,",");
      
      numsur = atoi(variables);
      printf("numsur = %i\n", numsur);
      
      variables = strtok(NULL," \t");
      variables = strtok(NULL," \t");
      variables = strtok(NULL,",");
      
      n = atoi(variables);
      printf("n = %i\n", n);
    }
  }

  //********************** MEMORY ALLOCATIONS ****************************************************
  
  dim = n;
  hlbrt1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim);
  data1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim);
  phase_node = (float **) malloc(numnodos * sizeof(float*));  
  for (i=0;i<numnodos;i++) phase_node[i] = (float *) malloc(dim * sizeof(float));
  phase_surr = (float ***) malloc(numnodos * sizeof(float**));    
  for (i=0;i<numnodos;i++){
    phase_surr[i] = (float **) malloc(numsur * sizeof(float *));
    for (j=0;j<numsur;j++) phase_surr[i][j] = (float *) malloc(dim * sizeof(float));
  }
  correlsurr = (float *) malloc(numsur * sizeof(float));
  dim2 = dim+1;
  correlsurrtime = (float *) malloc(dim2 * sizeof(float));
  correltime = (float *) malloc(dim2 * sizeof(float));
  correlnodes = (float **) malloc(numnodos * sizeof(float *));
  for (i=0;i<numnodos;i++) correlnodes[i] = (float *) malloc(numnodos * sizeof(float));
  delaynodes = (int **) malloc(numnodos * sizeof(int*));
  for (i=0;i<numnodos;i++) delaynodes[i] = (int *) malloc(numnodos * sizeof(int));

  pf = fftw_plan_dft_1d(dim, hlbrt1, hlbrt1, FFTW_FORWARD, FFTW_ESTIMATE);  // in place transforms: from hlbrt1 -> hlbrt1
  pb = fftw_plan_dft_1d(dim, hlbrt1, hlbrt1, FFTW_BACKWARD, FFTW_ESTIMATE); // in place transforms: from hlbrt1 -> hlbrt1

  nodo = 0;
  while (nodo < numnodos){
    cnt = 0; // "cnt" conta surrogates do nodo "nodo"

    //********************* READS CHANNEL NODO ********************************************************
  
    if(getline(&line, &line_size, fp) != -1 && line[0] != '#'){
      i=0;
      variables = strtok(line," \t");
      aux = atof(variables);
      data1[i] = aux + 0.*I;
      while (i < dim-1){ 
	// while (variables) {
	variables = strtok(NULL,"\t");
	i++;
	aux = atof(variables);
	data1[i] = aux + 0.*I;
      }
      fflush( stdout );

      // prepare for in place fftw, initialize hlbrt1
      for (i=0;i<dim;i++){
	hlbrt1[i] = data1[i];
      }

      // fftw mapped onto hlbrt1
      fftw_execute(pf);

      // set DC component and Nyquist component to zero
      hlbrt1[0] = 0.+I*0.;
      index = floor(dim/2);
      hlbrt1[index] = 0.+I*0.;
      
      // multiply positive frequency components by -j
      for (i=1;i<dim/2;i++){
	hlbrt1[i] = cimag(hlbrt1[i]) - I*creal(hlbrt1[i]);
      }
      
      // multiply negative frequency components by +j
      for (i=dim/2+1; i<dim; i++){
	hlbrt1[i] = -cimag(hlbrt1[i]) + I*creal(hlbrt1[i]);
      }

      // inverse fftw mapped onto hlbrt1
      fftw_execute(pb);

      // normalization of inverse fourier transform:: hlbrt1 is Hilbert transform of data1
      for (i=0;i<dim;i++) hlbrt1[i]=(float)(hlbrt1[i]/dim);

      // compute signal phase
      for (i=0;i<dim;i++){
	phase_node[nodo][i] = atan2(creal(hlbrt1[i]),creal(data1[i]));
      }
    }


    //********************* READS CHANNEL NODO's SURROGATES *******************************************
    
    while (getline(&line, &line_size, fp) != -1 && line[0] != '#' && cnt < numsur){
      cnt++;
      i=0;
      variables = strtok(line," \t");
      aux = atof(variables);
      data1[i] = aux + 0.*I;
      while (i < dim-1){ 
	// while (variables) {
	variables = strtok(NULL,"\t");
	i++;
	aux = atof(variables);
	data1[i] = aux + 0.*I;
      }
      fflush( stdout );

      // prepare for in place fftw, initialize hlbrt1
      for (i=0;i<dim;i++)   hlbrt1[i] = data1[i];
    
      // fftw mapped onto hlbrt1
      fftw_execute(pf);
  
      // set DC component and Nyquist component to zero
      hlbrt1[0] = 0.+I*0.;
      index = floor(dim/2);
      hlbrt1[index] = 0.+I*0.;
      
      //multiply positive frequency components by -j
      for (i=1;i<dim/2;i++)	hlbrt1[i] = cimag(hlbrt1[i]) - I*creal(hlbrt1[i]);
      
      //multiply negative frequency components by +j
      for (i=dim/2+1; i<dim; i++)	hlbrt1[i] = -cimag(hlbrt1[i]) + I*creal(hlbrt1[i]);
 
      // inverse fftw mapped onto hlbrt1
      fftw_execute(pb);

      // normalization of inverse fourier transform:: hlbrt1 is Hilbert transform of data1
      for (i=0;i<dim;i++)      hlbrt1[i]=(float)(hlbrt1[i]/dim);

      // compute signal phase for surrogate 'cnt-1'
      for (i=0;i<dim;i++){
	phase_surr[nodo][cnt-1][i] = atan2(creal(hlbrt1[i]),creal(data1[i]));
      }
    }
    nodo++;
  }
  fclose(fp);
  
  //Computation of phase correlations
  index = floor(dim/2);
  for (i=0;i<numnodos-1;i++){
    //    for (k=i+1;k<numnodos;k++){
    for (k=i+1;k<numnodos;k++){
      // computes phase correlation between node "i" and surrogate "cnt-1" of node "k" for different delays -dim/2 < delay < dim/2
      for (cnt=0;cnt<numsur;cnt++){
	printf("Correlation %i\t%i\t%i\n",i,k,cnt);
	for (delay=0;delay<=dim;delay++)    correlsurrtime[delay] = 0.+0.*I;
      
	// Use equation (14) Schmidt, Petkov, Richardson, Terry, "Dynamics on networks:..." for delay >= 0
	for (delay=0;delay<=index;delay++){
	  for (j=0;j<dim-delay;j++){
	    correlsurrtime[delay+index] += cos(phase_node[i][j+delay]-phase_surr[k][cnt][j])+I*sin(phase_node[i][j+delay]-phase_surr[k][cnt][j]);
	  }
	  correlsurrtime[delay+index]/=(float)(dim-delay);
	  correlsurrtime[delay+index] = sqrt(creal(correlsurrtime[delay+index])*creal(correlsurrtime[delay+index])+cimag(correlsurrtime[delay+index])*cimag(correlsurrtime[delay+index]));
	}
	// Use equation (14) Schmidt, Petkov, Richardson, Terry, "Dynamics on networks:..." for delay < 0
	for (delay=-index;delay<0;delay++){
	  for (j=abs(delay);j<dim;j++){
	    correlsurrtime[delay+index] += cos(phase_node[i][j+delay]-phase_surr[k][cnt][j])+I*sin(phase_node[i][j+delay]-phase_surr[k][cnt][j]);
	  }
	  correlsurrtime[delay+index]/=(float)(dim-abs(delay));
	  correlsurrtime[delay+index] = sqrt(creal(correlsurrtime[delay+index])*creal(correlsurrtime[delay+index])+cimag(correlsurrtime[delay+index])*cimag(correlsurrtime[delay+index]));
	}

	// maximize correlsurrtime over delay:: the maximum correlation will be assigned as correlsur
	correlsurr[cnt] = correlsurrtime[0];
	for (delay=1;delay<=dim;delay++){
	  if (correlsurrtime[delay] > correlsurr[cnt])    correlsurr[cnt] = correlsurrtime[delay];
	}
      }
      
      // Compute histogram of correlations between node "i" and surrogates of node "k"
      for (l=0;l<nbin;l++) hist[l] = 0.;
      maxx = correlsurr[0];
      minn = correlsurr[0];
      for (l=1;l<numsur;l++){
	if (correlsurr[l] > maxx)   maxx = correlsurr[l];
	if (correlsurr[l] < minn)   minn = correlsurr[l];
	}
      dx = (float)((maxx-minn)/nbin);
      for (l=0;l<numsur;l++){
	j = floor((correlsurr[l]-minn)/dx);
	hist[j] += (float)(1./numsur/dx);
      }
      
      // We compute threshold by computing the integral of hist until a bin just below 0.95, then do a linear regression to find intermediate point (correlatoin threshold) for which integral is 0.95
      prob = 0.;
      l=0;
      while (prob<0.95){ // 95% confidence level
	prob += hist[l];
	l++;
      }
      l--;
      if (prob>0.95){
	prob -= hist[l];
	// linear regression:
	threshold = (l-1)*dx+(0.95-prob)/(hist[l+1]-hist[l])*dx;
      }
      else threshold = l*dx;

      // computes phase correlation between node "k" and surrogate "cnt-1" of node "i" for different delays -dim/2 < delay < dim/2
      for (cnt=0;cnt<numsur;cnt++){
	printf("Correlation %i\t%i\t%i\n",k,i,cnt);
      	for (delay=0;delay<=dim;delay++)    correlsurrtime[delay] = 0.+0.*I;
	
	// Use equation (14) Schmidt, Petkov, Richardson, Terry, "Dynamics on networks:..." for delay >= 0
	for (delay=0;delay<=index;delay++){
	  for (j=0;j<dim-delay;j++){
	    correlsurrtime[delay+index] += cos(phase_node[k][j+delay]-phase_surr[i][cnt][j])+I*sin(phase_node[k][j+delay]-phase_surr[i][cnt][j]);
	  }
	  correlsurrtime[delay+index]/=(float)(dim-delay);
	  correlsurrtime[delay+index] = sqrt(creal(correlsurrtime[delay+index])*creal(correlsurrtime[delay+index])+cimag(correlsurrtime[delay+index])*cimag(correlsurrtime[delay+index]));
	}
	// Use equation (14) Schmidt, Petkov, Richardson, Terry, "Dynamics on networks:..." for delay < 0
	for (delay=-index;delay<0;delay++){
	  for (j=abs(delay);j<dim;j++){
	    correlsurrtime[delay+index] += cos(phase_node[k][j+delay]-phase_surr[i][cnt][j])+I*sin(phase_node[k][j+delay]-phase_surr[i][cnt][j]);
	  }
	  correlsurrtime[delay+index]/=(float)(dim-abs(delay));
	  correlsurrtime[delay+index] = sqrt(creal(correlsurrtime[delay+index])*creal(correlsurrtime[delay+index])+cimag(correlsurrtime[delay+index])*cimag(correlsurrtime[delay+index]));
	}
	
	// maximize correlsurrtime over delay:: the maximum correlation will be assigned as correlsur
	correlsurr[cnt] = correlsurrtime[0];
	for (delay=1;delay<=dim;delay++){
	  if (correlsurrtime[delay] > correlsurr[cnt])    correlsurr[cnt] = correlsurrtime[delay];
	}
      }
      
      // Compute histogram of correlations between node "i" and surrogates of node "k"
      for (l=0;l<nbin;l++) hist[l] = 0.;
      maxx = correlsurr[0];
      minn = correlsurr[0];
      for (l=1;l<numsur;l++){
	if (correlsurr[l] > maxx)   maxx = correlsurr[l];
	if (correlsurr[l] < minn)   minn = correlsurr[l];
      }
      dx = (float)((maxx-minn)/nbin);
      for (l=0;l<numsur;l++){
	j = floor((correlsurr[l]-minn)/dx);
	hist[j] += (float)(1./numsur/dx);
      }
      
      // We compute threshold by computing the integral of hist until a bin just below 0.95, then do a linear regression to find intermediate point (correlatoin threshold) for which integral is 0.95
      prob = 0.;
      l=0;
      while (prob<0.95){ // 95% confidence level
	prob += hist[l];
	l++;
      }
      l--;
      if (prob>0.95){
	prob -= hist[l];
	// linear regression:
	threshold2 = (l-1)*dx+(0.95-prob)/(hist[l+1]-hist[l])*dx;
      }
      else threshold2 = l*dx;
      
      if (threshold2>threshold) threshold=threshold2;
      
      //**************
      
      // compute cross-correlations between nodes from "phase_nodes" and compare with "threshold"s
      for (delay=0;delay<=dim;delay++)    correltime[delay] = 0.+0.*I;
      for (delay=0;delay<=index;delay++){
	for (j=0;j<dim-delay;j++){
	  correltime[delay+index] += cos(phase_node[i][j+delay]-phase_node[k][j])+I*sin(phase_node[i][j+delay]-phase_node[k][j]);
	}
	correltime[delay+index]/=(float)(dim-delay);
	correltime[delay+index] = sqrt(creal(correltime[delay+index])*creal(correltime[delay+index])+cimag(correltime[delay+index])*cimag(correltime[delay+index]));
	}
	for (delay=-index;delay<0;delay++){
	  for (j=abs(delay);j<dim;j++){
	    correltime[delay+index] += cos(phase_node[i][j+delay]-phase_node[k][j])+I*sin(phase_node[i][j+delay]-phase_node[k][j]);
	  }
	  correltime[delay+index]/=(float)(dim-abs(delay));
	  correltime[delay+index] = sqrt(creal(correltime[delay+index])*creal(correltime[delay+index])+cimag(correltime[delay+index])*cimag(correltime[delay+index]));
	}
	
	// maximize correltime over delay:: the maximum correlation will be assigned as correlnodes
	correlnodes[i][k] = correltime[0];
	delaynodes[i][k] = -index;
	for (delay=1;delay<=dim;delay++){
	  if (correltime[delay] > correlnodes[i][k]){    
	    correlnodes[i][k] = correltime[delay];
	    correlnodes[k][i] = correlnodes[i][k];
	    delaynodes[i][k] = delay-index;
	  }
	}
      
	// filter correlnodes according to threshold of correlation of node "i" and surrogates of "k"
	if (correlnodes[i][k] < threshold){
	  correlnodes[i][k] = 0.;
	  correlnodes[k][i] = 0.;
	}
    }
  }
    
  // filter correlations directionally (according to sign of delaynodes[i][k])
  for (i=0;i<numnodos-1;i++){
    for (k=i+1;k<numnodos;k++){
      if (correlnodes[i][k] != 0.){
  	if (delaynodes[i][k] == 0){
  	  correlnodes[i][k] = 0.;
  	  correlnodes[k][i] = 0.;
  	}
  	if (delaynodes[i][k] < 0){ // i influences k; k is influenced by i.
	  correlnodes[i][k] = 0.;
  	}
  	if (delaynodes[i][k] > 0){ // k influences i; i is influenced by k.
  	  correlnodes[k][i] = 0.;
  	}
      }
    }
  }

  // filter correlations according to three nodes comparison
  for (i=0;i<numnodos;i++){
    for (j=0;j<numnodos;j++){
      for (k=0;k<numnodos;k++){
  	if (correlnodes[i][k] > correlnodes[i][j] && correlnodes[k][j] > correlnodes[i][j]){
  	  correlnodes[i][j]=0.;
  	  break;
  	}
      }
    }
  }

  // filter correlations according to four nodes comparison
  if (four == 1){
    for (i=0;i<numnodos;i++){
      for (j=0;j<numnodos;j++){
	for (k=0;k<numnodos;k++){
	  if (correlnodes[i][k] > correlnodes[i][j]){
	    for (l=0;l<numnodos;l++){
	      if (correlnodes[k][l] > correlnodes[i][j] && correlnodes[l][j] > correlnodes[i][j]){
		correlnodes[i][j]=0.;
		break;
	      }
	    }
	    break;
	  }
	}
      }
    }
  }
  
  // print functional network matrix
  fp = fopen("functional_matrix.dat","w");
  for (i=0;i<numnodos;i++){
    for (j=0;j<numnodos;j++)      fprintf(fp,"%i\t%i\t%f\n",i,j,correlnodes[i][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  int G[numnodos][numnodos];
  
  // print associated pairs
  fp = fopen("association_nodes.dat","w");
  for (i=0;i<numnodos;i++){
    for (j=0;j<numnodos;j++)
      if(correlnodes[i][j] != 0){
	fprintf(fp,"%i\t%i\n",i,j);
	G[i][j] = 1;
      }
      else G[i][j] = 0;
  }
  fclose(fp);

  
  
  fftw_destroy_plan(pf);  fftw_destroy_plan(pb);
  fftw_free(data1); fftw_free(hlbrt1);
  
  return 0;
}
