//Compile linking with fftw3:
// gcc <NAMEOFFILE>.c -lfftw3 -lm

//Comparison between Discrete Hilbert Transform with Analytic Hilbert Transform (in gnuplot, for instance): p "<nameoffile>.dat" u 1:2 w l, "" u 1:3 w l

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

//*******************************************************************************************************************************************************


int main(){
  int i,j,k,T,n,cnt,nodo,numnodos,numsur,dim,index;
  float *s,*phi1,*phi2,**correltime,**correlnodes,*correlsurr,*correlsurrtime,*correllow,*correlhigh,aux,pi;
  FILE *fp;
  char *line,*variables;
  size_t line_size=1024;
  fftw_complex *in, *out,  *data1,*data2,*hlbrt1,*hlbrt2,*analytic_signal;
  fftw_plan pf,pb;

  pi = acos(-1.);
  
  line = malloc(line_size * sizeof(char));
  variables = malloc(line_size * sizeof(char));
  
  int fa = 500;

  
  //  dimensao de hlbrt1:: **Não** precisa ser potencia de 2
  dim = 100000;
  hlbrt1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim);

  pf = fftw_plan_dft_1d(dim, hlbrt1, hlbrt1, FFTW_FORWARD, FFTW_ESTIMATE); // in place fft transform
  pb = fftw_plan_dft_1d(dim, hlbrt1, hlbrt1, FFTW_BACKWARD, FFTW_ESTIMATE);

  //*********************************************************************************
  // ********************************  COSINE  **************************************
  //*********************************************************************************
  
  for (i=0;i<dim;i++){
    hlbrt1[i] = cos(2.*pi*(float)(1.*i/fa))+I*0.; // fa = freq amostragem
  }
    
  fftw_execute(pf); 
  
  fp = fopen("out_fftcos.dat","w");
  for(i=0;i<dim;i++) fprintf(fp,"%lf\t%lf\n",creal(hlbrt1[i]),cimag(hlbrt1[i]));
  fclose(fp);

  hlbrt1[0] = 0.+I*0.;
  index = floor(dim/2);
  hlbrt1[index] = 0.+I*0.;
  
  //multiply positive frequency components by -j
  for (i=1;i<dim/2;i++){
    hlbrt1[i] = cimag(hlbrt1[i])-I*creal(hlbrt1[i]);
  }

  //multiply negative frequency components by +j
  for (i=dim/2+1; i<dim; i++){
    hlbrt1[i] = -cimag(hlbrt1[i])+I*creal(hlbrt1[i]);
  }
 
  fftw_execute(pb);

  for (i=0;i<dim;i++) hlbrt1[i]=(float)(hlbrt1[i]/dim);
 
  fp = fopen("out_htcos.dat","w");
  for(i=0;i<dim;i++)  fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(float)(1.*i/fa),creal(hlbrt1[i]),sin(2.*pi*(float)(1.*i/fa)),cimag(hlbrt1[i]));
  fclose(fp);

  
  //*********************************************************************************
  // ********************************  SINE  ****************************************
  //*********************************************************************************
  
  for (i=0;i<dim;i++){
    hlbrt1[i] = sin(2.*pi*(float)(1.*i/fa))+I*0.; // fa = freq amostragem
  }
    
  fftw_execute(pf); 
  
  fp = fopen("out_fftsin.dat","w");
  for(i=0;i<dim;i++) fprintf(fp,"%lf\t%lf\n",creal(hlbrt1[i]),cimag(hlbrt1[i]));
  fclose(fp);

  hlbrt1[0] = 0.+I*0.;
  index = floor(dim/2);
  hlbrt1[index] = 0.+I*0.;
  
  //multiply positive frequency components by -j
  for (i=1;i<dim/2;i++){
    hlbrt1[i] = cimag(hlbrt1[i])-I*creal(hlbrt1[i]);
  }

  //multiply negative frequency components by +j
  for (i=dim/2+1; i<dim; i++){
    hlbrt1[i] = -cimag(hlbrt1[i])+I*creal(hlbrt1[i]);
  }
 
  fftw_execute(pb);

  for (i=0;i<dim;i++) hlbrt1[i]=(float)(hlbrt1[i]/dim);
 
  fp = fopen("out_htsin.dat","w");
  for(i=0;i<dim;i++)  fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(float)(1.*i/fa),creal(hlbrt1[i]),-cos(2.*pi*(float)(1.*i/fa)),cimag(hlbrt1[i]));
  fclose(fp);

  
  //*********************************************************************************
  // ********************************  SINC  ****************************************
  //*********************************************************************************

  hlbrt1[0]=0.+I*0.;
  for (i=1;i<dim;i++){
    hlbrt1[i] = sin((float)(1.*i/fa))/(float)(1.*i/fa)+I*0.;
  }
    
  fftw_execute(pf); /* repeat as needed */
  
  fp = fopen("out_fftsinc.dat","w");
  for(i=0;i<dim;i++) fprintf(fp,"%lf\t%lf\n",creal(hlbrt1[i]),cimag(hlbrt1[i]));
  fclose(fp);

  hlbrt1[0] = 0.+I*0.;
  index = floor(dim/2);
  hlbrt1[index] = 0.+I*0.;
  
  //multiply positive frequency components by -j
  for (i=1;i<dim/2;i++){
    hlbrt1[i] = cimag(hlbrt1[i])-I*creal(hlbrt1[i]);
  }
  
  //multiply negative frequency components by +j
  for (i=dim/2+1; i<dim; i++){
    hlbrt1[i] = -cimag(hlbrt1[i])+I*creal(hlbrt1[i]);
  }
 
  fftw_execute(pb);

  for (i=0;i<dim;i++) hlbrt1[i]=(float)(hlbrt1[i]/dim);
 
  fp = fopen("out_htsinc.dat","w");
  for(i=0;i<dim;i++)  fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(float)(1.*i/fa),creal(hlbrt1[i]),pi*sin(0.5*(float)(1.*i/fa))*sin(0.5*(float)(1.*i/fa))/0.5/(float)(1.*i/fa),cimag(hlbrt1[i]));
  fclose(fp);

  //*********************************************************************************
  //***************************  Function 1/(1+x^2)  ********************************
  //*********************************************************************************
  
  for (i=0;i<dim;i++){
    hlbrt1[i] = (float)(1./(1.+(float)(1.*i*i/fa/fa)))+I*0.; // fa = freq amostragem
  }
    
  fftw_execute(pf); 
  
  fp = fopen("out_fftcauchy.dat","w");
  for(i=0;i<dim;i++) fprintf(fp,"%lf\t%lf\n",creal(hlbrt1[i]),cimag(hlbrt1[i]));
  fclose(fp);

  hlbrt1[0] = 0.+I*0.;
  index = floor(dim/2);
  hlbrt1[index] = 0.+I*0.;
  
  //multiply positive frequency components by -j
  for (i=1;i<dim/2;i++){
    hlbrt1[i] = cimag(hlbrt1[i])-I*creal(hlbrt1[i]);
  }

  //multiply negative frequency components by +j
  for (i=dim/2+1; i<dim; i++){
    hlbrt1[i] = -cimag(hlbrt1[i])+I*creal(hlbrt1[i]);
  }
 
  fftw_execute(pb);

  for (i=0;i<dim;i++) hlbrt1[i]=(float)(hlbrt1[i]/dim);
 
  fp = fopen("out_htcauchy.dat","w");
  for(i=0;i<dim;i++)      fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(float)(1.*i/fa),creal(hlbrt1[i]),(float)(1.*i/fa)/(1.+(float)(1.*i/fa)*(float)(1.*i/fa)),cimag(hlbrt1[i]));
  fclose(fp);
    
  fftw_destroy_plan(pf);  fftw_destroy_plan(pb);
  fftw_free(hlbrt1); 
  
  return 0;
  
}
