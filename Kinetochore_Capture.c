*****************************************************************************************************************************************************************************
CKAP5 regulates microtubule-chromosome attachment fidelity by stabilizing CENP-E at kinetochores CKAP5 stabilizes CENP-E at kinetochores by regulating microtubule-chromosome attachment
*****************************************************************************************************************************************************************************
*  School of Mathematical and Computational Sciences, Indian Association for the Cultivation of Science, Jadavpur, Kolkata 700032, India
*  School of Biology, Indian Institute of Science Education and Research, Thiruvananthapuram, Vithura, Thiruvananthapuram 695551, Kerala, India
*  Department of Human Molecular Genetics and Biochemistry, Faculty of Medicine, Tel Aviv University, Tel Aviv, Israel
*  New Chemistry Unit, Jawaharlal Nehru Centre for Advanced Scientific Research, Bengaluru 560064, Karnataka, India.
*  (c) 2024
*
*  DESIGN & DEVELOPMENTS
*
*  Pinaki Nayak, Apurba Sarkar, Raja Paul
*  School of Mathematical and Computational Sciences, Indian Association for the Cultivation of Science, Jadavpur, Kolkata 700032, India
*
*
*  R. Bhagya Lakshmi, Vishnu M. Nair, Delvin K. Pauly, S. Sajana, Sanusha M. G., Tapas K Manna
*  School of Biology, Indian Institute of Science Education and Research, Thiruvananthapuram, Vithura, Thiruvananthapuram 695551, Kerala, India
*
*  Linoy Raz, Uri Ben-David
*  Department of Human Molecular Genetics and Biochemistry, Faculty of Medicine, Tel Aviv University, Tel Aviv, Israel
*
*  Akshay Saroha, Pratibha Kumari, Sarit S. Agasti
*  New Chemistry Unit, Jawaharlal Nehru Centre for Advanced Scientific Research, Bengaluru 560064, Karnataka, India.
*
*
*
*****************************************************************************************************************************************************************************

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


//..........................................parameters................................................//

double pi = 3.14159265359;	//value of pi

//cell size, major and minor axes (um)//
double A_x;// = 10.0;
double A_y;// = 10.0;
double A_z;// = 10.0;

long idum;	//random number seed
int idum2;

double R_ch = 0.5;	//CH-arm radius um
double R_kt = 0.5;	//Radius of KT (um)
double R_cMT_visual = 0.1;	//cMT cylinder radius for povray visualization (um)

double R_cent = 0.5;	//centrosome radius (um)

double originx = 0.0;	//cell center
double originy = 0.0;
double originz = 0.0;

int N_cent;	//number of centrosomes
int N_chr;	//number of chromosomes
int Num_cMT;	//number of cMT per centrosome, cMT==>centrosomal MT
int Max_capture; // maximum MTs one KT can capture

#define N_cent 2
#define N_chr  80
#define Num_cMT 500
#define Max_capture 5
//dynamic instability parameters//

double fc = 0.025340;	//catastrophe frequency (1/sec), ~3.0*vg/(4*Rcell)
double fr = 0.003;	//rescue frequency (1/sec)
double fc_l = 0.08;
double fc_stall = 0.05;
double fs = 1.7;
double vg = 0.25;	//growth velocity (um/sec)
double vs = 0.4;	//shrinkage velocity (um/sec)
double mean = 0.0;
double standard_dev  = 0.3;
double D_kt = 0.0001;
double mu = 5.0;

/////////////////////////////////
double KT_MT_int_rad = 2*R_kt;
double kt_kt_l0 = 2*(R_ch + R_kt);
double kt_mt_l0 = R_kt;
double kt_kt_ksp = 100.0;
double kt_mt_ksp = 10.0;
double kt_fibrils_sp = 2.5;
double P_motor = 0.5;
double A_pe = 25.0;
double L_pe = 3.0;
double A_steric = 0.1;
double A_boundary = 10.0;
double L_boundary = 3.0;
double detach_motor_distance = 25.0;
double k_syn = 0.001;
double k_mero = 0.001;
double k1_mero = 1.0;
double k2_mero = 1.0;

int t;	//iteration timestep starting from zero
int t_tot = 240000;	//total number of simulation timestep; (t_tot*h/60) min
double h = 0.01;	//delta t in sec
double angle_h = 0.01;
int iteration_upper_limit = 5000000;	//maximum number of trials to place chromosomes
int move_upper_limit = 5000;
int detach_relax = 1;
int ens_counter; //iteratively simulating ensemble 1 to N_ens
int N_ens; //total number of ensembles
int N_start;

//.....................................memory allocation.....................................//

///////////////////////////////////////
char fName[50];
FILE *fp;

char fName1[50];
FILE *fp1;

char gName[50];
FILE *gp;

double cent_x[N_cent],cent_y[N_cent],cent_z[N_cent];	//centrosome coordinates

double chrom_x[N_chr],chrom_y[N_chr],chrom_z[N_chr];	//chromosome center coordinates
double chrom_theta[N_chr],chrom_phi[N_chr];

double KTl_x[N_chr],KTl_y[N_chr],KTl_z[N_chr];	//left KT coordinates
double KTr_x[N_chr],KTr_y[N_chr],KTr_z[N_chr];	//right KT coordinates
double label_KTl[N_chr],label_KTr[N_chr];	//label = 1.0 ==> unattached KT, -1.0 ==> KT attached to MT
int capture_KTl_cent[N_chr][Max_capture],capture_KTr_cent[N_chr][Max_capture]; //KT captured by centrosome
int capture_KTl_MT[N_chr][Max_capture],capture_KTr_MT[N_chr][Max_capture]; //KT captured by MT

double force_KTl[N_chr][3],force_KTr[N_chr][3],force_chr[N_chr][3]; //force on the KTs

double l_cmt[N_cent][Num_cMT],cMTtheta[N_cent][Num_cMT],cMTphi[N_cent][Num_cMT],label_cMT[N_cent][Num_cMT];	//cMT length, direction, label_cMT= +1.0==>unattached MT,=-1.0==>attached to KT
double cMTtipx[N_cent][Num_cMT],cMTtipy[N_cent][Num_cMT],cMTtipz[N_cent][Num_cMT];	//cMT tip coordinates
int capture_MT_chrom[N_cent][Num_cMT],capture_MT_KT[N_cent][Num_cMT];

double chrom_x_old,chrom_y_old,chrom_z_old;
double KTl_x_old,KTl_y_old,KTl_z_old,KTr_x_old,KTr_y_old,KTr_z_old;

double cMTflag[N_cent][Num_cMT]; //cMTflag = 1.0 ==>growing cMT, -1.0 ==>shrinking cMT

double P_cat; //1-exp(-fc*h)
double P_res; //1-exp(-fr*h)


//..................................random number generator ran2()...................................//


void nrerror0(const char error_text[])
{
        printf("Numerical Recipes run-time error...\n");
        printf("%s\n",error_text);
        printf("...now exiting to system...\n");
        exit(1);
}

/**********************************************************/

/***** random number generator ran2 ***********************/
#include <math.h>

#define M 714025
#define IA 1366
#define IC 150889

double ran2(int *idum)

{
    static long iy,irn[98];
    static int iff=0;
    int j;
    void nrerror();

    if (*idum < 0 || iff == 0) {
        iff=1;
        if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
        for (j=1;j<=97;j++) {
            *idum=(IA*(*idum)+IC) % M;
            irn[j]=(*idum);
        }
        *idum=(IA*(*idum)+IC) % M;
        iy=(*idum);
    }
    j=1 +(int)( 97.0*iy/M);
    if (j > 97 || j < 1) nrerror0("RAN2: This cannot happen.");
    iy=irn[j];
    *idum=(IA*(*idum)+IC) % M;
    irn[j]=(*idum);
    return (double)iy/M;
}

#undef M
#undef IA
#undef IC
/****************************************************************/



//...............................gaussian random number generator.............................//

#define NITER 4

void psdes(unsigned long *lword, unsigned long *irword)
{
unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
static unsigned long c1[NITER]={
0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
static unsigned long c2[NITER]={
0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
for (i=0;i<NITER;i++) {
  ia=(iswap=(*irword)) ^ c1[i];
itmpl = ia & 0xffff;
itmph = ia >> 16;
ib=itmpl*itmpl+ ~(itmph*itmph);
*irword=(*lword) ^ (((ia = (ib >> 16) |
((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
 *lword=iswap;
 }
}

float ran4(long *idum)
{
void psdes(unsigned long *lword, unsigned long *irword);
unsigned long irword,itemp,lword;
static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
static unsigned long jflone = 0x00004080;
static unsigned long jflmsk = 0xffff007f;
#else
static unsigned long jflone = 0x3f800000;
static unsigned long jflmsk = 0x007fffff;
#endif
if (*idum < 0) {
idums = -(*idum);
*idum=1;
}
irword=(*idum);
lword=idums;
psdes(&lword,&irword);
itemp=jflone | (jflmsk & irword);
++(*idum);
return (*(float *)&itemp)-1.0;
}

float gasdev(long *idum,double mean,double standard_deviation)
{
  float ran4(long *idum);
  static int iset=0;
  static double y1,y2;
  float fac,rsq,v1,v2;
  if (*idum < 0) iset=0;
  if (iset == 0) {
    do {
      v1=2.0*ran4(idum)-1.0;
      v2=2.0*ran4(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    y1=v1*fac;
    y2=v2*fac;
    iset=1;
    return (mean+standard_deviation*(double)y2);
  } else {
    iset=0;
    return (mean+standard_deviation*(double)y1);
  }
}

//.........................................place components.................................//

///////////////////////////////////////////////////////
double distance(double x1,double y1,double z1,double x2,double y2,double z2){
       double d;
       d=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
       return d;
}
////////////////////////////////////////////////////////
double place_centrosome(){
       int i;

       for(i=0;i<N_cent;i++){
          cent_x[i] = -2.0*i*0.6*A_x+0.6*A_x;
	  			cent_y[i] = 0.0;
	  			cent_z[i] = 0.0;
          //printf("%f %f %f\n",cent_x[i],cent_y[i],cent_z[i]);
       }
       //printf("centrosomes placed");
    }
//////////////////////////////////////////////////
double norm(double x, double y, double z){
  return(sqrt(x*x+y*y+z*z));
}
/////////////////////////////////////////////////
int OutOfSpheroid(double x,double y,double z,double a, double b, double c){
  //Here, we assume the spheroid is centered at 0,0,0 and
  //a,b,c are the radii in 3 directions

  if(((x*x)/(a*a)+(y*y)/(b*b)+(z*z)/(c*c))<1.0){
    return(0);//inside the spheroid
  }else return(1);//outside the spheroid
}
///////////////////////////////////////////////////
int OutOfCell(int chr_count){

    int chr_tag,KTl_tag,KTr_tag;
    chr_tag = OutOfSpheroid(chrom_x[chr_count],chrom_y[chr_count],chrom_z[chr_count],A_x,A_y,A_z);
    KTl_tag = OutOfSpheroid(KTl_x[chr_count],KTl_y[chr_count],KTl_z[chr_count],A_x,A_y,A_z);
    KTr_tag = OutOfSpheroid(KTr_x[chr_count],KTr_y[chr_count],KTr_z[chr_count],A_x,A_y,A_z);

    if(chr_tag == 1 || KTl_tag == 1 || KTr_tag == 1){
      return 1;
    }
    else{
      return 0;
    }
}
//////////////////////////////////////////////////
double get_chr_distance(int i,int j){
       double d,dx,dy,dz;

       dx = (chrom_x[i] - chrom_x[j]);
       dy = (chrom_y[i] - chrom_y[j]);
       dz = (chrom_z[i] - chrom_z[j]);

       d = norm(dx,dy,dz);

       return(d);
}
///////////////////////////////////////////////////
double avoid_chr_overlap(int chr_index){
       int i;
       double d,overlap_flag;
       overlap_flag=0.0;
       for(i=0;i<chr_index;i++){

         //check chr overalap
         d=distance(chrom_x[i],chrom_y[i],chrom_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
         if(d<2*(R_ch)){
           overlap_flag=1.0;
         }
         d=distance(KTl_x[i],KTl_y[i],KTl_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
         if(d<(R_ch+R_kt)){
           overlap_flag=1.0;
         }
         d=distance(KTr_x[i],KTr_y[i],KTr_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
         if(d<(R_ch+R_kt)){
           overlap_flag=1.0;
         }
         //check left kt overlap
         d=distance(KTl_x[i],KTl_y[i],KTl_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
         if(d<2*(R_kt)){
           overlap_flag=1.0;
         }
         d=distance(KTr_x[i],KTr_y[i],KTr_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
         if(d<2*(R_kt)){
           overlap_flag=1.0;
         }
         d=distance(chrom_x[i],chrom_y[i],chrom_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
         if(d<(R_ch+R_kt)){
           overlap_flag=1.0;
         }
         //check right kt overlap
         d=distance(KTl_x[i],KTl_y[i],KTl_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
         if(d<2*(R_kt)){
           overlap_flag=1.0;
         }
         d=distance(KTr_x[i],KTr_y[i],KTr_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
         if(d<2*(R_kt)){
           overlap_flag=1.0;
         }
         d=distance(chrom_x[i],chrom_y[i],chrom_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
         if(d<(R_ch+R_kt)){
           overlap_flag=1.0;
         }

       }
       return overlap_flag;
}
/////////////////////////////////////////////////////////////////
double rotate_chr_kt(int chr_index,double theta, double phi){
       int i;
       double r;

       r = R_ch+R_kt;

       KTl_x[chr_index] = chrom_x[chr_index] + r*sin(theta)*cos(phi);
       KTl_y[chr_index] = chrom_y[chr_index] + r*sin(theta)*sin(phi);
       KTl_z[chr_index] = chrom_z[chr_index] + r*cos(theta);

       KTr_x[chr_index] = chrom_x[chr_index] - r*sin(theta)*cos(phi);
       KTr_y[chr_index] = chrom_y[chr_index] - r*sin(theta)*sin(phi);
       KTr_z[chr_index] = chrom_z[chr_index] - r*cos(theta);

}
/////////////////////////////////////////////////////////////////////
double place_chromosomes(){
       int i,chr_count,attempt_count,out_of_sph_tag;
       double x,y,z,theta,phi,overlap_tag;

       chr_count = 0;
       attempt_count = 0;

       while(chr_count<N_chr){

         chrom_x[chr_count] = 2*A_x*ran2(& idum2) - A_x;
         chrom_y[chr_count] = 2*A_y*ran2(& idum2) - A_y;
         chrom_z[chr_count] = 2*A_z*ran2(& idum2) - A_z;

         theta = pi*ran2(& idum2);
         phi = 2*pi*ran2(& idum2);
         chrom_theta[chr_count] = theta;
         chrom_phi[chr_count] = phi;

         rotate_chr_kt(chr_count,theta,phi);
         overlap_tag = avoid_chr_overlap(chr_count);
         out_of_sph_tag = OutOfCell(chr_count);

         if(overlap_tag != 1.0 && out_of_sph_tag != 1){
           chr_count++;
           attempt_count = 0;
         }

         attempt_count++;
         if(attempt_count == iteration_upper_limit){
           printf("Exceeded maximum trials to place chromosomes, exiting\n");
           exit(0);
         }

       }
       //printf("chromosomes placed");


}
/////////////////////////////////////////////////////
double save_old_position(int i){

			 chrom_x_old = chrom_x[i];
			 chrom_y_old = chrom_y[i];
			 chrom_z_old = chrom_z[i];

			 KTl_x_old = KTl_x[i];
			 KTl_y_old = KTl_y[i];
			 KTl_z_old = KTl_z[i];

			 KTr_x_old = KTr_x[i];
			 KTr_y_old = KTr_y[i];
			 KTr_z_old = KTr_z[i];

}

//////////////////////////////////////////////////////

double copy_old_position(int i){

			 chrom_x[i]= chrom_x_old;
		   chrom_y[i]= chrom_y_old;
	     chrom_z[i]= chrom_z_old;

			 KTl_x[i] = KTl_x_old;
	     KTl_y[i] = KTl_y_old;
	     KTl_z[i] = KTl_z_old;

	     KTr_x[i] = KTr_x_old;
	     KTr_y[i] = KTr_y_old;
	     KTr_z[i] = KTr_z_old;
}

//////////////////////////////////////////////////////////

double avoid_chr_overlap_move(int chr_index){
  int i;
  double d,overlap_flag;
  overlap_flag=0.0;
  for(i=0;i<N_chr;i++){

    if(i==chr_index){
      continue;
    }
    //check chr overalap
    d=distance(chrom_x[i],chrom_y[i],chrom_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
    if(d<2*(R_ch)){
      overlap_flag=1.0;
    }
    d=distance(KTl_x[i],KTl_y[i],KTl_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
    if(d<(R_ch+R_kt)){
      overlap_flag=1.0;
    }
    d=distance(KTr_x[i],KTr_y[i],KTr_z[i],chrom_x[chr_index],chrom_y[chr_index],chrom_z[chr_index]);
    if(d<(R_ch+R_kt)){
      overlap_flag=1.0;
    }
    //check left kt overlap
    d=distance(KTl_x[i],KTl_y[i],KTl_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
    if(d<2*(R_kt)){
      overlap_flag=1.0;
    }
    d=distance(KTr_x[i],KTr_y[i],KTr_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
    if(d<2*(R_kt)){
      overlap_flag=1.0;
    }
    d=distance(chrom_x[i],chrom_y[i],chrom_z[i],KTl_x[chr_index],KTl_y[chr_index],KTl_z[chr_index]);
    if(d<(R_ch+R_kt)){
      overlap_flag=1.0;
    }
    //check right kt overlap
    d=distance(KTl_x[i],KTl_y[i],KTl_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
    if(d<2*(R_kt)){
      overlap_flag=1.0;
    }
    d=distance(KTr_x[i],KTr_y[i],KTr_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
    if(d<2*(R_kt)){
      overlap_flag=1.0;
    }
    d=distance(chrom_x[i],chrom_y[i],chrom_z[i],KTr_x[chr_index],KTr_y[chr_index],KTr_z[chr_index]);
    if(d<(R_ch+R_kt)){
      overlap_flag=1.0;
    }

  }
  return overlap_flag;
}
/////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

double move_chromosomes(){
       int i,chr_count,attempt_count,out_of_sph_tag;
       double d,dx,dy,dz,theta,phi,overlap_tag;

       chr_count = 0;
       attempt_count = 0;

       while(chr_count<N_chr){
         save_old_position(chr_count);

         dx = sqrt(2*D_kt)*gasdev(& idum,mean,standard_dev);
         dy = sqrt(2*D_kt)*gasdev(& idum,mean,standard_dev);
         dz = sqrt(2*D_kt)*gasdev(& idum,mean,standard_dev);

         //update coordianates by forces

         KTl_x[chr_count] += (1/mu)*force_KTl[chr_count][0]*h + dx;
         KTl_y[chr_count] += (1/mu)*force_KTl[chr_count][1]*h + dy;
         KTl_z[chr_count] += (1/mu)*force_KTl[chr_count][2]*h + dz;

         KTr_x[chr_count] += (1/mu)*force_KTr[chr_count][0]*h + dx;
         KTr_y[chr_count] += (1/mu)*force_KTr[chr_count][1]*h + dy;
         KTr_z[chr_count] += (1/mu)*force_KTr[chr_count][2]*h + dz;

         chrom_x[chr_count] =  (KTl_x[chr_count] + KTr_x[chr_count])/2.0;
         chrom_y[chr_count] =  (KTl_y[chr_count] + KTr_y[chr_count])/2.0;
         chrom_z[chr_count] =  (KTl_z[chr_count] + KTr_z[chr_count])/2.0;

         overlap_tag = avoid_chr_overlap_move(chr_count);
         out_of_sph_tag = OutOfCell(chr_count);

        if(overlap_tag == 1.0 || out_of_sph_tag == 1){
          copy_old_position(chr_count);

        }

        // update coordinates by diffusion

        KTl_x[chr_count] +=  dx;
        KTl_y[chr_count] +=  dy;
        KTl_z[chr_count] +=  dz;

        KTr_x[chr_count] +=  dx;
        KTr_y[chr_count] +=  dy;
        KTr_z[chr_count] +=  dz;

        chrom_x[chr_count] += dx;
        chrom_y[chr_count] += dy;
        chrom_z[chr_count] += dz;

        overlap_tag = avoid_chr_overlap_move(chr_count);
        out_of_sph_tag = OutOfCell(chr_count);

       if(overlap_tag == 1.0 || out_of_sph_tag == 1){
         copy_old_position(chr_count);

       }
        chr_count++;


    }
}
///////////////////////////////////////////////////////////

//..........................initialize arrays....................................//

/////////////////////////////////////
double cMT_direction(int i){
      int j,out_of_sph_tag;
      double x,y,z;
      double theta,phi;
      for(j=0;j<Num_cMT;j++){
         if(label_cMT[i][j]==1.0){
           while(1<2){
                theta = pi*ran2(& idum2);
                phi = 2.0*pi*ran2(& idum2);
                x = cent_x[i]+0.01*sin(theta)*cos(phi);
                y = cent_y[i]+0.01*sin(theta)*sin(phi);
                z = cent_z[i]+0.01*cos(theta);
                out_of_sph_tag = OutOfSpheroid(x,y,z,A_x,A_y,A_z);
                if(out_of_sph_tag==0){
                  break;
                }
           }
           cMTtheta[i][j] = theta;
           cMTphi[i][j]  = phi;
           l_cmt[i][j] = 0.0;
         }
      }
   }
///////////////////////////////////
////////////////////////////////
double init_MT_length(){
      int i,j;
      for(i=0;i<N_cent;i++){
         for(j=0;j<Num_cMT;j++){
            l_cmt[i][j] = 0.0;
         }
      }
   }
////////////////////////////////
/////////////////////////////////////
double init_flag(){
      int i,j,k;
      for(i=0;i<N_cent;i++){
         for(j=0;j<Num_cMT;j++){
            cMTflag[i][j] = 1.0;  //growing cMT
            label_cMT[i][j] = 1.0; //free cMT, not attached to any KT
         }
      }
      for(k=0;k<N_chr;k++){
         label_KTl[k] = 1.0;
         label_KTr[k] = 1.0;
      }
   }
////////////////////////////////////////

////////////////////////////////////////
double Initialize_Arrays(){
       int i,j;
       init_MT_length();
       init_flag();
       for(i=0;i<N_cent;i++){
         cMT_direction(i);
         for(j=0;j<Num_cMT;j++){
           capture_MT_KT[i][j]=-1;
           capture_MT_chrom[i][j]=-1;
         }
      }
      for(i=0;i<N_chr;i++){
        label_KTl[i]=1.0;
        label_KTr[i]=1.0;
        for(j=0;j<Max_capture;j++){

            capture_KTl_MT[i][j]=-1;
            capture_KTr_MT[i][j]=-1;
            capture_KTl_cent[i][j]=-1;
            capture_KTr_cent[i][j]=-1;

        }

      }
    }
/////////////////////////////////////////


//............................................microtubule dynamics.............................//

//////////////////////////////////////////
double initiate_newMT(int i,int j){
       int out_of_sph_tag;
       double x,y,z;
       double theta,phi;
       while(1<2){

            theta = pi*ran2(& idum2);
            phi = 2.0*pi*ran2(& idum2) - pi;

            x = cent_x[i]+0.01*sin(theta)*cos(phi);
            y = cent_y[i]+0.01*sin(theta)*sin(phi);
            z = cent_z[i]+0.01*cos(theta);
            out_of_sph_tag = OutOfSpheroid(x,y,z,A_x,A_y,A_z);
            if(out_of_sph_tag==0){
              cMTtheta[i][j] = theta;
              cMTphi[i][j]  = phi;
              l_cmt[i][j] = 0.0;
              cMTflag[i][j] =-1.0;
              break;
            }
       }
       //printf("new mt %d	%d %d\n",t,i,j);
}
////////////////////////////////////////

///////////////////////////////////////////
double mt_depoly_force(int i,int j){
  double d,d1,d2,f;
  int chr,kt;
  chr = capture_MT_chrom[i][j];
  kt = capture_MT_KT[i][j];
  if(kt==0){
      d = distance(cMTtipx[i][j],cMTtipy[i][j],cMTtipz[i][j],KTl_x[chr],KTl_y[chr],KTl_z[chr]);
    //  d1 = distance(KTl_x[chr],KTl_y[chr],KTl_z[chr],cent_x[i],cent_y[i],cent_z[i]);
  }
  if(kt==1){
      d = distance(cMTtipx[i][j],cMTtipy[i][j],cMTtipz[i][j],KTr_x[chr],KTr_y[chr],KTr_z[chr]);
      //d1 = distance(KTr_x[chr],KTr_y[chr],KTr_z[chr],cent_x[i],cent_y[i],cent_z[i]);
  }

  f = -kt_mt_ksp*(d - kt_mt_l0);
  //d2 = distance(cMTtipx[i][j],cMTtipy[i][j],cMTtipz[i][j],cent_x[i],cent_y[i],cent_z[i]);

  if(f<0.0){
      return abs(f);
  }
  else{
    return 0.0;
  }

}
////////////////////////////////////////////////
double cMTshrink(int i,int j){
      double p,q,f_res,f_depoly;
      q = 0;
      //printf("%d\n",flag);
      if(cMTflag[i][j]==-1.0){
        if(l_cmt[i][j]>(vs*h)){ //restricitng MT length to be positive
          l_cmt[i][j]=l_cmt[i][j]-vs*h;
        }else{q=1;}
      }
      p = ran2(& idum2);

      if(label_cMT[i][j] == -1.0){

        f_depoly = mt_depoly_force(i,j);
        if(f_depoly>0.0){
          f_res = 1.0 - exp(-f_depoly/fs);
          if(f_res<fr){
            f_res = fr;
          }

        }
        else{
          f_res = fr;
        }
        P_res = 1.0 - exp(-f_res*h);
      }
      else{
        P_res = 1.0-exp(-fr*h);
      }

      if(q==1){
        initiate_newMT(i,j);
        p=-1.0;
      }
      if(cMTflag[i][j]==-1.0 && p<=P_res){
        cMTflag[i][j] = 1.0;
      }
   }
////////////////////////////////////////

//////////////////////////////////////
double cMTgrow(int i,int j){
      double p,f_c,f_load;
      if(cMTflag[i][j]==1.0){
        l_cmt[i][j] = l_cmt[i][j]+vg*h;
        p = ran2(& idum2);
        if(label_cMT[i][j] == -1.0){
          f_c = fc_l*l_cmt[i][j];
          //f_load = mt_depoly_force(i,j);
        //  f_c=fc_stall/(1+(fc_stall/fc -1)*exp(-f_load/fs));
          P_cat = 1.0-exp(-f_c*h);
        }
        else{
          P_cat = 1.0-exp(-fc*h);
        }
        if(p<P_cat){
          cMTflag[i][j] = -1.0;
        } 	//state switching
      }
   }
///////////////////////////////////////

double Grow_Shrink_cMT(){
       int i,j;
       for(i=0;i<N_cent;i++){
          for(j=0;j<Num_cMT;j++){

               cMTshrink(i,j);
               cMTgrow(i,j);

          }
       }
    }
//////////////////////////////////////
double cMT_tip(){
      int i,j,k;
      for(i=0;i<N_cent;i++){
         for(j=0;j<Num_cMT;j++){
            cMTtipx[i][j] = cent_x[i] + l_cmt[i][j]*sin(cMTtheta[i][j])*cos(cMTphi[i][j]);
            cMTtipy[i][j] = cent_y[i] + l_cmt[i][j]*sin(cMTtheta[i][j])*sin(cMTphi[i][j]);
            cMTtipz[i][j] = cent_z[i] + l_cmt[i][j]*cos(cMTtheta[i][j]);

         }
      }
   }
/////////////////////////////////////
double rotate_MT(){
    int i,j;
    double x,y,z,x0,y0,z0,dx,dy,dz,l,theta_kt,phi_kt;
    int ch,kt;
    for(i=0;i<N_cent;i++){
      for(j=0;j<Num_cMT;j++){
        if(label_cMT[i][j] == -1.0){
          ch = capture_MT_chrom[i][j];
          kt = capture_MT_KT[i][j];
          if(kt == 0){
              x = KTl_x[ch];
              y = KTl_y[ch];
              z = KTl_z[ch];
          }
          if(kt == 1){
              x = KTr_x[ch];
              y = KTr_y[ch];
              z = KTr_z[ch];
          }

          x0 = cent_x[i];
          y0 = cent_y[i];
          z0 = cent_z[i];

          dx = x - x0;
          dy = y - y0;
          dz = z - z0;

          l= norm(dx,dy,dz);
          if(l<0.001){
            l=0.001;
          }

          if(abs(dx)<0.001){
            if(dx>0.0){
                dx =0.001;
            }
            else{
              dx =-0.001;
            }

          }

          if(abs(dy)<0.001){
            if(dy>0.0){
                dy =0.001;
            }
            else{
              dy =-0.001;
            }

          }
          if(abs(dz)<0.001){
            if(dz>0.0){
                dz =0.001;
            }
            else{
              dz =-0.001;
            }

          }

          phi_kt = atan2(dy,dx);
          theta_kt = acos(dz/l);

          cMTtheta[i][j] += angle_h*(theta_kt - cMTtheta[i][j]);
          cMTphi[i][j] += angle_h*(phi_kt - cMTphi[i][j]);

        }
      }
    }
}
//////////////////////////////////
double Check_cMT_boundary_encounter(){
       int i,j,out_of_sph_tag;
       double x,y,z;
       for(i=0;i<N_cent;i++){
          for(j=0;j<Num_cMT;j++){
             if(cMTflag[i][j]==1.0){
               x = cMTtipx[i][j];
               y = cMTtipy[i][j];
               z = cMTtipz[i][j];
               out_of_sph_tag = OutOfSpheroid(x,y,z,A_x,A_y,A_z);
               if(out_of_sph_tag==1){
                 cMTflag[i][j] = -1.0; //MT catastrophe upon hitting boundary
               }
             }
          }
       }
   }
/////////////////////////////////////

//.............................................kinetochore capture subroutines.....................................//
///////////////////////////////////////////////
double Shortest_Distance(double x1, double y1, double z1, double x2, double y2, double z2, double x, double y, double z){

       //This routine returns the shortest distance from the point x,y,z
       //on the line joining x1,y1,z1 to x2,y2,z2

       //cross product components of P2P1 x PP1
       double P2P1xP1P_1 = (y2-y1)*(z-z1)-(z2-z1)*(y-y1);
       double P2P1xP1P_2 = (z2-z1)*(x-x1)-(x2-x1)*(z-z1);
       double P2P1xP1P_3 = (x2-x1)*(y-y1)-(y2-y1)*(x-x1);

       //modulus of the cross product
       double P2P1xP1P = sqrt(P2P1xP1P_1*P2P1xP1P_1+P2P1xP1P_2*P2P1xP1P_2+P2P1xP1P_3*P2P1xP1P_3);
       //modulus of P2P1
       double P2P1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
       double d = P2P1xP1P/P2P1;

       return(d);
}
///////////////////////////////////////////////
int check_capture_present(int kt,int chr,int cent,int mt,int max){
    int i,flag;
    flag =0;
    for(i=0;i<Max_capture;i++){
      if(kt==0){
        if((capture_KTl_cent[chr][i] == cent) && (capture_KTl_MT[chr][i] == mt)){
          flag=1;
          break;
        }
      }
      if(kt==1){
        if((capture_KTr_cent[chr][i] == cent) && (capture_KTr_MT[chr][i] == mt)){
          flag=1;
          break;
        }
      }
    }
    return flag;
}

//////////////////////////////////////////////////
double check_amphitelic(int chr){
    int i,j;
    double flag;
    int KTl_cs1,KTl_cs2,KTr_cs1,KTr_cs2;
    KTl_cs1 = 0;
    KTl_cs2 = 0;
    KTr_cs1 = 0;
    KTr_cs2 = 0;

    for(i=0;i<Max_capture;i++){
      if(capture_KTl_cent[chr][i]==0){
        KTl_cs1++;
      }
      if(capture_KTl_cent[chr][i]==1){
        KTl_cs2++;
      }
      if(capture_KTr_cent[chr][i]==0){
        KTr_cs1++;
      }
      if(capture_KTr_cent[chr][i]==1){
        KTr_cs2++;
      }
    }

    //unattached
    if(KTl_cs1 == 0 && KTl_cs2 == 0 && KTr_cs1 == 0 && KTr_cs2 == 0){
      flag = 3.0; //unattached
      return flag;
    }
    //merotelic check
    if(KTl_cs1 > 0 && KTl_cs2 >0){
      flag = 1.0; //merotelic
      return flag;
    }
    if(KTr_cs1 > 0 && KTr_cs2 >0){
      flag = 1.0; //merotelic
      return flag;
    }

    //syntelic check
    if(KTl_cs1 > 0 && KTl_cs2 == 0 && KTr_cs1 > 0 && KTr_cs2 == 0){
      flag = 2.0; //syntelic
      return flag;
    }
    if(KTl_cs1 == 0 && KTl_cs2 > 0 && KTr_cs1 == 0 && KTr_cs2 > 0){
      flag = 2.0; //syntelic
      return flag;
    }

    //amphitelic check
    if(KTl_cs1 > 0 && KTl_cs2 == 0 && KTr_cs1 == 0 && KTr_cs2 > 0){
      flag = 0.0; //amphitelic
      return flag;
    }
    if(KTl_cs1 == 0 && KTl_cs2 > 0 && KTr_cs1 > 0 && KTr_cs2 == 0){
      flag = 0.0; //amphitelic
      return flag;
    }

    return 3.0;

}
////////////////////////////////////////////////////
double check_MT_KT_encounter(){
       int i,j,k,cap_count,cent,mt,flag_capture_present,flag_amphitelic;
       double d,p;
       double kt_cap_cent_old[Max_capture],kt_cap_mt_old[Max_capture];

       for(i=0;i<N_chr;i++){

         flag_amphitelic = check_amphitelic(i);
         if(flag_amphitelic == 0.0){
          // printf("amphitelic!");
           continue;
         }

         cap_count = 0;

         for(j=0;j<Max_capture;j++){
           kt_cap_mt_old[j] = -1;
           kt_cap_cent_old[j] = -1;
         }


         if(label_KTl[i] == -1.0){
           for(j=0;j<Max_capture;j++){
             cent = capture_KTl_cent[i][j];
             mt = capture_KTl_MT[i][j];
             if(cent!=-1 && mt!=-1){
               d = distance(KTl_x[i],KTl_y[i],KTl_z[i],cMTtipx[cent][mt],cMTtipy[cent][mt],cMTtipz[cent][mt]);
               if(d < detach_motor_distance){
                 kt_cap_cent_old[cap_count] = cent;
                 kt_cap_mt_old[cap_count] = mt;
                 cap_count++;
                 capture_MT_chrom[cent][mt] = i;
                 capture_MT_KT[cent][mt] = 0;
               }
               else{
                 label_cMT[cent][mt] = 1.0;
               }
             }
           }
         }

         for(j=0;j<Max_capture;j++){
           capture_KTl_cent[i][j] = kt_cap_cent_old[j];
           capture_KTl_MT[i][j] = kt_cap_mt_old[j];
         }

         if(cap_count<Max_capture){
           for(j=0;j<N_cent;j++){
             for(k=0;k<Num_cMT;k++){
               flag_capture_present = check_capture_present(0,i,j,k,cap_count);

               if(flag_capture_present == 0 && cap_count < Max_capture && label_cMT[j][k] == 1.0){
                 d = distance(KTl_x[i],KTl_y[i],KTl_z[i],cMTtipx[j][k],cMTtipy[j][k],cMTtipz[j][k]);
                 p = ran2(& idum2);
                 if(d < KT_MT_int_rad && p < P_motor){
                   label_KTl[i] = -1.0;
                   label_cMT[j][k] = -1.0;
                   capture_KTl_cent[i][cap_count] = j;
                   capture_KTl_MT[i][cap_count] = k;
                   cap_count++;
                   capture_MT_chrom[j][k] = i;
                   capture_MT_KT[j][k] = 0;
                 }
               }
             }
           }
         }
         if(cap_count == 0){
           label_KTl[i] = 1.0;
         }


         // now the same for right KT

         cap_count = 0;

         for(j=0;j<Max_capture;j++){
           kt_cap_mt_old[j] = -1;
           kt_cap_cent_old[j] = -1;
         }

         if(label_KTr[i] == -1.0){
           for(j=0;j<Max_capture;j++){
             cent = capture_KTr_cent[i][j];
             mt = capture_KTr_MT[i][j];
             if(cent!=-1 && mt!=-1){
               d = distance(KTr_x[i],KTr_y[i],KTr_z[i],cMTtipx[cent][mt],cMTtipy[cent][mt],cMTtipz[cent][mt]);
               if(d < detach_motor_distance){
                 kt_cap_cent_old[cap_count] = cent;
                 kt_cap_mt_old[cap_count] = mt;
                 cap_count++;
                 capture_MT_chrom[cent][mt] = i;
                 capture_MT_KT[cent][mt] = 1;
               }
               else{
                 label_cMT[cent][mt] = 1.0;
               }
             }
           }
         }

         for(j=0;j<Max_capture;j++){
           capture_KTr_cent[i][j] = kt_cap_cent_old[j];
           capture_KTr_MT[i][j] = kt_cap_mt_old[j];
         }

         if(cap_count<Max_capture){
           for(j=0;j<N_cent;j++){
             for(k=0;k<Num_cMT;k++){
               flag_capture_present = check_capture_present(1,i,j,k,cap_count);

               if(flag_capture_present == 0 && cap_count < Max_capture && label_cMT[j][k] == 1.0){
                 d = distance(KTr_x[i],KTr_y[i],KTr_z[i],cMTtipx[j][k],cMTtipy[j][k],cMTtipz[j][k]);
                 p = ran2(& idum2);
                 if(d < KT_MT_int_rad && p < P_motor){
                   label_KTr[i] = -1.0;
                   label_cMT[j][k] = -1.0;
                   capture_KTr_cent[i][cap_count] = j;
                   capture_KTr_MT[i][cap_count] = k;
                   cap_count++;
                   capture_MT_chrom[j][k] = i;
                   capture_MT_KT[j][k] = 1;
                 }
               }
             }
           }
         }
         if(cap_count == 0){
           label_KTr[i] = 1.0;
         }

       }
}
/////////////////////////////////////////////
void captured_ch_print(int ens,int t){
  int i;
  FILE *fp;
  char name[50];

  sprintf(name,"ensemble_%d_capture_count_at_%d.dat",ens,t);
  fp = fopen(name,"w");
  int c_amp,c_syn,c_mero;
  c_amp = 0;
  c_syn = 0;
  c_mero = 0;

  for(i=0;i<N_chr;i++){
    if(check_amphitelic(i)==0.0){
      c_amp++;
    }
    if(check_amphitelic(i)==1.0){
      c_mero++;
    }
    if(check_amphitelic(i)==2.0){
      c_syn++;
    }


  }
  fprintf(fp,"%d %d %d\n",c_amp,c_syn,c_mero);
  fclose(fp);

}

void capture_mt_count(int ens,int t){
  FILE *fp;
  char name[50];
  int KTl_cs1,KTl_cs2,KTr_cs1,KTr_cs2,chr,i;
  sprintf(name,"ensemble_%d_MT_capture_count_%d.dat",ens,t);
  fp = fopen(name,"w");

  for(chr=0;chr<N_chr;chr++){
    KTl_cs1 = 0;
    KTl_cs2 = 0;
    KTr_cs1 = 0;
    KTr_cs2 = 0;
    for(i=0;i<Max_capture;i++){
      if(capture_KTl_cent[chr][i]==0){
        KTl_cs1++;
      }
      if(capture_KTl_cent[chr][i]==1){
        KTl_cs2++;
      }
      if(capture_KTr_cent[chr][i]==0){
        KTr_cs1++;
      }
      if(capture_KTr_cent[chr][i]==1){
        KTr_cs2++;
      }
    }
    fprintf(fp,"%d %d %d %d %d\n",chr, KTl_cs1,KTl_cs2,KTr_cs1,KTr_cs2);
  }
	fclose(fp);
}
/////////////////////////////////////////////////////

//....................................forces.................................//

double KT_MT_forces(){
       int i,j,cent,mt;
       double d,fx,fy,fz,dx,dy,dz,kt_mt_force;

       for(i=0;i<N_chr;i++){

         fx = 0.0;
         fy = 0.0;
         fz = 0.0;

         if(label_KTl[i] == -1.0){                       // force on left KT due to attached MTs
           for(j=0;j<Max_capture;j++){
             cent =  capture_KTl_cent[i][j];
             mt = capture_KTl_MT[i][j];

             if(cent != -1 && mt != -1 ){

               dx = KTl_x[i] - cMTtipx[cent][mt];
               dy = KTl_y[i] - cMTtipy[cent][mt];
               dz = KTl_z[i] - cMTtipz[cent][mt];
               d = norm(dx,dy,dz);
               if(d<0.001){
                 d=0.001;
               }
               if(d > R_kt){
                 kt_mt_force = -1.0*kt_mt_ksp*(d-kt_mt_l0);
                }
                else{
                  kt_mt_force = -1.0*kt_mt_ksp*(d-kt_mt_l0);
                }

               fx += kt_mt_force*dx/d;
               fy += kt_mt_force*dy/d;
               fz += kt_mt_force*dz/d;

             }
           }
           force_KTl[i][0] += fx;
           force_KTl[i][1] += fy;
           force_KTl[i][2] += fz;
         }



         fx = 0.0;
         fy = 0.0;
         fz = 0.0;

         if(label_KTr[i] == -1.0){                     // force on right KT due to attached MTs
           for(j=0;j<Max_capture;j++){
             cent =  capture_KTr_cent[i][j];
             mt = capture_KTr_MT[i][j];

             if(cent != -1 && mt != -1){

               dx = KTr_x[i] - cMTtipx[cent][mt];
               dy = KTr_y[i] - cMTtipy[cent][mt];
               dz = KTr_z[i] - cMTtipz[cent][mt];
               d = norm(dx,dy,dz);

               if(d<0.001){
                 d=0.001;
               }
               if(d > R_kt){
                 kt_mt_force = -1.0*kt_mt_ksp*(d-kt_mt_l0);
                }
                else{
                  kt_mt_force = -1.0*kt_mt_ksp*(d-kt_mt_l0);
                }
               fx += kt_mt_force*dx/d;
               fy += kt_mt_force*dy/d;
               fz += kt_mt_force*dz/d;

             }
           }
           force_KTr[i][0] += fx;
           force_KTr[i][1] += fy;
           force_KTr[i][2] += fz;
         }



       }
   }
/////////////////////////////////////////////////////////////////////
double KT_fibrils_forces(){
       int i,j,k;
       double d,dx,dy,dz,f,fx,fy,fz;

       for(i=0;i<N_chr;i++){
         for(j=0;j<N_cent;j++){
           for(k=0;k<Num_cMT;k++){

             dx = KTl_x[i] - cMTtipx[j][k];
             dy = KTl_y[i] - cMTtipy[j][k];
             dz = KTl_z[i] - cMTtipz[j][k];
             d = norm(dx,dy,dz);

             if(d < R_kt){
               f = kt_fibrils_sp*(R_kt-d);

               dx = KTl_x[i] - cent_x[j];
               dy = KTl_y[i] - cent_y[j];
               dz = KTl_z[i] - cent_z[j];
               d = norm(dx,dy,dz);

               fx = f*dx/d;
               fy = f*dy/d;
               fz = f*dz/d;

               force_KTl[i][0] += fx;
               force_KTl[i][1] += fy;
               force_KTl[i][2] += fz;
             }

             dx = KTr_x[i] - cMTtipx[j][k];
             dy = KTr_y[i] - cMTtipy[j][k];
             dz = KTr_z[i] - cMTtipz[j][k];
             d = norm(dx,dy,dz);

             if(d < R_kt){
               f = kt_fibrils_sp*(R_kt-d);

               dx = KTl_x[i] - cent_x[j];
               dy = KTl_y[i] - cent_y[j];
               dz = KTl_z[i] - cent_z[j];
               d = norm(dx,dy,dz);

               fx = f*dx/d;
               fy = f*dy/d;
               fz = f*dz/d;

               force_KTr[i][0] += fx;
               force_KTr[i][1] += fy;
               force_KTr[i][2] += fz;
             }

           }
         }
       }

}
////////////////////////////////////////////////
double KT_KT_forces(){
       int i;
       double d,dx,dy,dz,fx,fy,fz,f;

       for(i=0;i<N_chr;i++){
            // forces due to stretching of sister KTs

           dx = KTl_x[i] - KTr_x[i];
           dy = KTl_y[i] - KTr_y[i];
           dz = KTl_z[i] - KTr_z[i];
           d = norm(dx,dy,dz);

           f = -kt_kt_ksp*(d - kt_kt_l0);
           fx = f*dx/d;
           fy = f*dy/d;
           fz = f*dz/d;

           force_KTl[i][0] += fx;
           force_KTl[i][1] += fy;
           force_KTl[i][2] += fz;

           force_KTr[i][0] -= fx;
           force_KTr[i][1] -= fy;
           force_KTr[i][2] -= fz;

       }
}
/////////////////////////////////////////////////////////////////
double chr_chr_steric_forces(){
       int i,j;
       double d,dx,dy,dz,f,fx,fy,fz;
       for(i=0;i<N_chr;i++){
         for(j=0;j<N_chr;j++){
           if(i==j){
             continue;
           }
           else{
             dx = chrom_x[i] - chrom_x[j];
             dy = chrom_y[i] - chrom_y[j];
             dz = chrom_z[i] - chrom_z[j];
             d = norm(dx,dy,dz);

             f = A_steric/(d*d);
             fx = f*dx/d;
             fy = f*dy/d;
             fz = f*dz/d;

             force_KTl[i][0] += fx;
             force_KTl[i][1] += fy;
             force_KTl[i][2] += fz;

             force_KTr[i][0] += fx;
             force_KTr[i][1] += fy;
             force_KTr[i][2] += fz;

           }
         }
       }


}
//////////////////////////////////////////////////
double distance_from_boundary(double x,double y,double z){
       double d,dx,dy,dz,cos_x,cos_y,cos_z,d_boundary,a;
       double x1,y1,z1;
       dx = x - originx;
       dy = y - originy;
       dz = z - originz;

       d = norm(dx,dy,dz);
       cos_x = dx/d;
       cos_y = dy/d;
       cos_z = dz/d;

       a = 1.0/sqrt((cos_x*cos_x)/(A_x*A_x) + (cos_y*cos_y)/(A_y*A_y) + (cos_z*cos_z)/(A_z*A_z));
       x1 = a*x;
       y1 = a*y;
       z1 = a*z;

       return distance(x,y,z,x1,y1,z1);

}
//////////////////////////////////////////////////
double chr_boundary_steric_force(){
       int i;
       double d,d_boundary,dx,dy,dz,f,fx,fy,fz;

       for(i=0;i<N_chr;i++){
         dx = chrom_x[i] - originx;
         dy = chrom_y[i] - originy;
         dz = chrom_z[i] - originz;

         d = norm(dx,dy,dz);
         d_boundary = distance_from_boundary(chrom_x[i],chrom_y[i],chrom_z[i]);

         f = A_boundary*exp(-1.0*d_boundary/L_boundary);

         fx = -f*dx/d;
         fy = -f*dy/d;
         fz = -f*dz/d;

         force_KTl[i][0] += fx;
         force_KTl[i][1] += fy;
         force_KTl[i][2] += fz;

         force_KTr[i][0] += fx;
         force_KTr[i][1] += fy;
         force_KTr[i][2] += fz;

       }
}
/////////////////////////////////////////////////
double chr_forces(){
       int i;
       double d;

       for(i=0;i<N_chr;i++){

         force_chr[i][0] += force_KTl[i][0] + force_KTr[i][0];
         force_chr[i][1] += force_KTl[i][1] + force_KTr[i][1];
         force_chr[i][2] += force_KTl[i][2] + force_KTr[i][2];

       }
}
double polar_ejection_force(){
       int i,j,k,flag;
       double d,dx,dy,dz,f,fx,fy,fz;

       for(i=0;i<N_chr;i++){

         for(j=0;j<N_cent;j++){
           flag = 0;

           for(k=0;k<Num_cMT;k++){
             d = distance(chrom_x[i],chrom_y[i],chrom_z[i],cMTtipx[j][k],cMTtipy[j][k],cMTtipz[j][k]);
             if(d<R_ch){
               flag = 1;
               break;
             }
           }

           if(flag == 1){

             dx = chrom_x[i] - cent_x[j];
             dy = chrom_y[i] - cent_y[j];
             dz = chrom_z[i] - cent_z[j];
             d = norm(dx,dy,dz);


             f = A_pe*exp(-d/L_pe);
             fx = f*dx/d;
             fy = f*dy/d;
             fz = f*dz/d;

             force_KTl[i][0] += fx;
             force_KTl[i][1] += fy;
             force_KTl[i][2] += fz;

             force_KTr[i][0] += fx;
             force_KTr[i][1] += fy;
             force_KTr[i][2] += fz;
           }
         }

      }
}
////////////////////////////////////////////////////////////////
double forces(){
    int i,j;
    for(i=0;i<N_chr;i++){
      for(j=0;j<3;j++){
        force_KTl[i][j] = 0.0;
        force_KTr[i][j] = 0.0;
        force_chr[i][j] = 0.0;
      }
    }

    KT_MT_forces();

    //KT_fibrils_forces();
    //chr_forces();
    KT_KT_forces();

    chr_chr_steric_forces();
    chr_boundary_steric_force();
    polar_ejection_force();
}

void print_forces(int t){
  int i,j;
  FILE *fp;
  char name[50];

  sprintf(name,"forces_at_%d.dat",t);
  fp=fopen(name,"w");
  for(i=0;i<N_chr;i++){
    fprintf(fp,"%lf %lf %lf\n",force_chr[i][0],force_chr[i][1],force_chr[i][2]);
  }
  fclose(fp);
}
/////////////////////////////////////////////////////////////////////////

//................................kinetochore-microtubule detachment dynamics.................................//
double correct_syntelic(int chr){
   int i;
   int cent,mt;
   double p,p_detach,syntelic_detach_rate,f_cs;
   for(i=0;i<Max_capture;i++){
     cent = capture_KTl_cent[chr][i];
     mt = capture_KTl_MT[chr][i];
     if(cent != -1){
       p = ran2(& idum2);
       f_cs = -kt_mt_ksp*(norm((cMTtipx[cent][mt]-KTl_x[cent]),(cMTtipy[cent][mt]-KTl_y[cent]),(cMTtipz[cent][mt]-KTl_z[cent])) - kt_mt_l0);
       syntelic_detach_rate = k_syn*exp(k1_mero*abs(f_cs));
       p_detach = 1.0 - exp(-syntelic_detach_rate*h);
       if(p<p_detach){
         capture_KTl_cent[chr][i] = -1;
         capture_KTl_MT[chr][i] = -1;
         label_cMT[cent][mt] = 1.0;
       }
     }
     cent = capture_KTr_cent[chr][i];
     mt = capture_KTr_MT[chr][i];
     if(cent != -1){
       p = ran2(& idum2);
       f_cs = -kt_mt_ksp*(norm((cMTtipx[cent][mt]-KTr_x[cent]),(cMTtipy[cent][mt]-KTr_y[cent]),(cMTtipz[cent][mt]-KTr_z[cent])) - kt_mt_l0);
       syntelic_detach_rate = k_syn*exp(k1_mero*abs(f_cs));
       p_detach = 1.0 - exp(-syntelic_detach_rate*h);
       if(p<p_detach){
         capture_KTr_cent[chr][i] = -1;
         capture_KTr_MT[chr][i] = -1;
         label_cMT[cent][mt] = 1.0;
       }
     }
   }
}
//////////////////////////////////////////////

double correct_merotelic(int chr){
  int i,KTl_cent,KTr_cent,cent,mt;
  double flag,p,p_detach,beta,cs1,cs2,n_cs,f_cs,fcs_x,fcs_y,fcs_z,f_kt,fkt_x,fkt_y,fkt_z,merotelic_detach_rate;

  // check left KT

  flag = 0.0;
  cs1=0.0;
  cs2=0.0;
  for(i=0;i<Max_capture;i++){
    if(capture_KTl_cent[chr][i]==0){
      cs1+=1.0;
    }
    if(capture_KTl_cent[chr][i]==1){
      cs2+=1.0;
    }
  }

  if(cs1>0 && cs2 >0){
    flag = 1.0;
  }

  if(flag == 1.0){
    for(i=0;i<Max_capture;i++){
      cent = capture_KTl_cent[chr][i];
      mt = capture_KTl_MT[chr][i];

      if(cent == -1){
        continue;
      }

      if(cent==0){
        n_cs=cs1/cs2;
      }
      else{
        n_cs=cs2/cs1;
      }

      f_cs = -kt_mt_ksp*(norm((cMTtipx[cent][mt]-KTl_x[cent]),(cMTtipy[cent][mt]-KTl_y[cent]),(cMTtipz[cent][mt]-KTl_z[cent])) - kt_mt_l0);
      fcs_x = f_cs*(KTl_x[cent]-cMTtipx[cent][mt])/norm((cMTtipx[cent][mt]-KTl_x[cent]),(cMTtipy[cent][mt]-KTl_y[cent]),(cMTtipz[cent][mt]-KTl_z[cent])) ;
      fcs_y = f_cs*(KTl_y[cent]-cMTtipy[cent][mt])/norm((cMTtipx[cent][mt]-KTl_x[cent]),(cMTtipy[cent][mt]-KTl_y[cent]),(cMTtipz[cent][mt]-KTl_z[cent])) ;
      fcs_z = f_cs*(KTl_z[cent]-cMTtipz[cent][mt])/norm((cMTtipx[cent][mt]-KTl_x[cent]),(cMTtipy[cent][mt]-KTl_y[cent]),(cMTtipz[cent][mt]-KTl_z[cent])) ;

      f_kt = -kt_kt_ksp*(norm(KTl_x[cent]-KTr_x[cent],KTl_y[cent]-KTr_y[cent],KTl_z[cent]-KTr_z[cent])-kt_kt_l0);
      fkt_x = f_kt*(KTr_x[cent]-KTl_x[cent])/norm(KTl_x[cent]-KTr_x[cent],KTl_y[cent]-KTr_y[cent],KTl_z[cent]-KTr_z[cent]);
      fkt_y = f_kt*(KTr_y[cent]-KTl_y[cent])/norm(KTl_x[cent]-KTr_x[cent],KTl_y[cent]-KTr_y[cent],KTl_z[cent]-KTr_z[cent]);
      fkt_z = f_kt*(KTr_z[cent]-KTl_z[cent])/norm(KTl_x[cent]-KTr_x[cent],KTl_y[cent]-KTr_y[cent],KTl_z[cent]-KTr_z[cent]);

      beta = (fkt_x*fcs_x + fkt_y*fcs_y + fkt_z*fcs_z)/(abs(f_cs)*abs(f_kt));

      merotelic_detach_rate = k_mero*n_cs*exp(k1_mero*abs(f_cs))*exp(-k2_mero*beta);
      p_detach = 1.0 - exp(-merotelic_detach_rate*h);

      p = ran2(& idum2);

      if(p<p_detach){
        capture_KTl_cent[chr][i] = -1;
        capture_KTl_MT[chr][i] = -1;
        label_cMT[cent][mt] = 1.0;
      }
    }
  }


  // check right KT

  flag = 0.0;
  cs1=0.0;
  cs2=0.0;
  for(i=0;i<Max_capture;i++){
    if(capture_KTr_cent[chr][i]==0){
      cs1+=1.0;
    }
    if(capture_KTr_cent[chr][i]==1){
      cs2+=1.0;
    }
  }

  if( cs1 > 0 && cs2 > 0 ){
    flag = 1.0;
  }

  if(flag == 1.0){
    for(i=0;i<Max_capture;i++){
      cent = capture_KTr_cent[chr][i];
      mt = capture_KTr_MT[chr][i];

      if(cent == -1){
        continue;
      }
      if(cent==0){
        n_cs=cs1/cs2;
      }
      else{
        n_cs=cs2/cs1;
      }

      f_cs = -kt_mt_ksp*(norm((cMTtipx[cent][mt]-KTr_x[cent]),(cMTtipy[cent][mt]-KTr_y[cent]),(cMTtipz[cent][mt]-KTr_z[cent])) - kt_mt_l0);
      fcs_x = f_cs*(KTl_x[cent]-cMTtipx[cent][mt])/norm((cMTtipx[cent][mt]-KTr_x[cent]),(cMTtipy[cent][mt]-KTr_y[cent]),(cMTtipz[cent][mt]-KTr_z[cent])) ;
      fcs_y = f_cs*(KTr_y[cent]-cMTtipy[cent][mt])/norm((cMTtipx[cent][mt]-KTr_x[cent]),(cMTtipy[cent][mt]-KTr_y[cent]),(cMTtipz[cent][mt]-KTr_z[cent])) ;
      fcs_z = f_cs*(KTr_z[cent]-cMTtipz[cent][mt])/norm((cMTtipx[cent][mt]-KTr_x[cent]),(cMTtipy[cent][mt]-KTr_y[cent]),(cMTtipz[cent][mt]-KTr_z[cent])) ;

      f_kt = -kt_kt_ksp*(norm(KTr_x[cent]-KTl_x[cent],KTr_y[cent]-KTl_y[cent],KTr_z[cent]-KTl_z[cent])-kt_kt_l0);
      fkt_x = f_kt*(KTl_x[cent]-KTr_x[cent])/norm(KTr_x[cent]-KTl_x[cent],KTr_y[cent]-KTl_y[cent],KTr_z[cent]-KTl_z[cent]);
      fkt_y = f_kt*(KTl_y[cent]-KTr_y[cent])/norm(KTr_x[cent]-KTl_x[cent],KTr_y[cent]-KTl_y[cent],KTr_z[cent]-KTl_z[cent]);
      fkt_z = f_kt*(KTl_z[cent]-KTr_z[cent])/norm(KTr_x[cent]-KTl_x[cent],KTr_y[cent]-KTl_y[cent],KTr_z[cent]-KTl_z[cent]);

      beta = (fkt_x*fcs_x + fkt_y*fcs_y + fkt_z*fcs_z)/(abs(f_cs)*abs(f_kt));

      merotelic_detach_rate = k_mero*n_cs*exp(k1_mero*abs(f_cs))*exp(-k2_mero*beta);
      p_detach = 1.0 - exp(-merotelic_detach_rate*h);

      p = ran2(& idum2);

      if(p<p_detach){
        capture_KTr_cent[chr][i] = -1;
        capture_KTr_MT[chr][i] = -1;
        label_cMT[cent][mt] = 1.0;
      }
    }
  }

}
//////////////////////////////
double check_capture_remains(int chr){
    int i,count_l,count_r;
    count_l = 0;
    count_r = 0;
    for(i=0;i<Max_capture;i++){
      if(capture_KTl_cent[chr][i]!=-1){
        count_l++;
      }
      if(capture_KTr_cent[chr][i]!=-1){
        count_r++;
      }
    }
    if(count_l == 0){
      label_KTl[chr] = 1.0;
    }
    if(count_r == 0){
      label_KTr[chr] = 1.0;
    }
}
//////////////////////////////////////
double correct_attachments(){
    int i,count;
    double flag;
    count = 0;
    while(count < detach_relax){

      for(i=0;i<N_chr;i++){
        if(label_KTl[i] == -1.0 || label_KTr[i] == -1.0){
          flag = check_amphitelic(i);
          if(flag == 1.0){
            correct_merotelic(i);
          }
          else if(flag == 2.0){
            correct_syntelic(i);
          }
          else{
            continue;
          }
        }
        check_capture_remains(i);
      }
      count++;
    }
}
////////////////////////////////////////////////
double print_chr_kt_coordinates(int ens, int t){
       int i;
       FILE *fp;
       char name[50];
       double d;
       sprintf(name,"Ensemble_%d_Chr_coordiantes_at_%d.dat",ens,t);
       fp = fopen(name,"w");
       //fp1 = fopen("kt_kt_distance.dat","a");
       for(i=0;i<N_chr;i++){
         d = norm(KTl_x[i]-KTr_x[i],KTl_y[i]-KTr_y[i],KTl_z[i]-KTr_z[i]);
         fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",check_amphitelic(i),chrom_x[i],chrom_y[i],chrom_z[i],KTl_x[i],KTl_y[i],KTl_z[i],KTr_x[i],KTr_y[i],KTr_z[i],d);

         //fprintf(fp1,"%lf\t",d);
       }

       //fprintf(fp1,"\n");
       //fclose(fp1);
       fclose(fp);

}
/////////////////////////////////////////////////////


//................................povray visual subroutines...................................//

/////////////////////////////////////////////
double init_povray_visual(){
      fprintf(fp,"#include \"colors.inc\"\n\n\n");
      fprintf(fp,"#include \"shapes.inc\"\n\n\n");
      fprintf(fp,"//Place the camera\n");
      fprintf(fp,"camera {\n");
      fprintf(fp,"sky <0,0,1>           //Don't change this\n");
      fprintf(fp,"direction <-1,0,0>    //Don't change this\n");
      fprintf(fp,"right <-4/3,0,0>      //Don't change this\n");
      fprintf(fp,"location <-0.0,-25.0,0.0> //Camera location\n");
      fprintf(fp,"look_at <0.0,0.0,0.0>     //Where camera is pointing\n");
      fprintf(fp,"angle 0      //Angle of the view--increase to see more, decrease to see less\n");
      fprintf(fp,"}\n");

      fprintf(fp,"//Ambient light to \"brighten up\" darker pictures\n");
      fprintf(fp,"global_settings { ambient_light White }\n");
      fprintf(fp,"//Set a background color\n");
      fprintf(fp,"background { color Black }\n");

      fprintf(fp,"light_source {\n");
      fprintf(fp," <10.0,0.0,0.0> \n");
      fprintf(fp,"color White*1.0  \n");
      fprintf(fp,"}\n");
}
////////////////////////////////////////
double Centrosome_povray_visual(double x1,double y1,double z1){
      fprintf(fp,"sphere{ <%f,%f,%f>,%f\n",x1,y1,z1,R_cent);
      fprintf(fp,"texture{ pigment{color rgbt <0.90, 0.0, 0.0, 0.0> }");
      fprintf(fp,"finish { phong 1 }\n");
      fprintf(fp,"}\n");

      fprintf(fp,"finish { ambient 0.5 specular 0.5 }\n");
      fprintf(fp,"no_shadow\n");
      fprintf(fp,"}\n"); //--------------------------------------------------------
}
////////////////////////////////////////
double Chromosome_arm_povray_visual(){
//chromosome
    int i;
    //fprintf(fp,"union\n");
    //fprintf(fp,"{\n");
    for (i=0;i<N_chr;i++){
      fprintf(fp,"object{cylinder{<%f,%f,%f>, <%f,%f,%f>, %f}\n",
	      KTl_x[i],KTl_y[i],KTl_z[i],KTr_x[i],KTr_y[i],KTr_z[i],0.2); //temp

      //fprintf(fp,"object{sphere{<%f,%f,%f>, %f}\n",
	     // chrom_x[i],chrom_y[i],chrom_z[i],R_ch);
        if(check_amphitelic(i)==0.0){
          fprintf(fp,"\n");
          fprintf(fp,"pigment{color rgb <1.0, 1.0, 0.0>}\n");
          fprintf(fp,"normal {bumps 1.0 scale 0.1}\n");
          fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
          fprintf(fp,"no_shadow\n");
          fprintf(fp,"\n");
          fprintf(fp,"}\n");
          fprintf(fp,"\n");
        }
        else if(check_amphitelic(i)==1.0){
          fprintf(fp,"\n");
          fprintf(fp,"pigment{color rgb <0.0, 1.0, 0.0>}\n");
          fprintf(fp,"normal {bumps 1.0 scale 0.1}\n");
          fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
          fprintf(fp,"no_shadow\n");
          fprintf(fp,"\n");
          fprintf(fp,"}\n");
          fprintf(fp,"\n");
        }
        else{
          fprintf(fp,"\n");
          fprintf(fp,"pigment{color rgb <1.0, 0.0, 1.0>}\n");
          fprintf(fp,"normal {bumps 1.0 scale 0.1}\n");
          fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
          fprintf(fp,"no_shadow\n");
          fprintf(fp,"\n");
          fprintf(fp,"}\n");
          fprintf(fp,"\n");
        }
     // fprintf(fp,"sphere{<%f,%f,%f>, %f}\n",
	     // chrom_x_tip2[i],chrom_y_tip2[i],chrom_z_tip2[i],R_ch);
    };
    /*fprintf(fp,"\n");
    fprintf(fp,"pigment{color rgb <1.0, 1.0, 0.0>}\n");
    fprintf(fp,"normal {bumps 1.0 scale 0.1}\n");
    fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
    fprintf(fp,"no_shadow\n");
    fprintf(fp,"\n");
    fprintf(fp,"}\n");
    fprintf(fp,"\n");*/
}
////////////////////////////////////////
double Kinetochore_povray_visual(){
      int i;
    //  fprintf(fp,"union\n");
    //  fprintf(fp,"{\n");
      for (i=0;i<N_chr;i++){
        //  fprintf(fp,"cylinder{<%f,%f,%f>, <%f,%f,%f>, %f}\n",
	      //KTl_x[i],KTl_y[i],KTl_z[i],chrom_x[i],chrom_y[i],chrom_z[i],R_kt);

          //fprintf(fp,"cylinder{<%f,%f,%f>, <%f,%f,%f>, %f}\n",
	      //KTr_x[i],KTr_y[i],KTr_z[i],chrom_x[i],chrom_y[i],chrom_z[i],R_kt);
          //fprintf(fp, "object{");
          fprintf(fp,"object{sphere{<%f,%f,%f>, %f}\n",KTl_x[i],KTl_y[i],KTl_z[i],R_kt);
           if(label_KTl[i] == -1){
             fprintf(fp,"pigment{color rgb <0.0, 0.0, 1.0>}\n");
             fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
             fprintf(fp,"no_shadow\n");
             fprintf(fp,"}\n");
             fprintf(fp,"\n");

           }
           else{
             fprintf(fp,"pigment{color rgb <0.0, 1.0, 1.0>}\n");
             fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
             fprintf(fp,"no_shadow\n");
              fprintf(fp,"}\n");
             fprintf(fp,"\n");
           }

          fprintf(fp,"object{sphere{<%f,%f,%f>, %f}\n",KTr_x[i],KTr_y[i],KTr_z[i],R_kt);
          if(label_KTr[i] == -1){
            fprintf(fp,"pigment{color rgb <0.0, 0.0, 1.0>}\n");
            fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
            fprintf(fp,"no_shadow\n");
            fprintf(fp,"}\n");
            fprintf(fp,"\n");

          }
          else{
            fprintf(fp,"pigment{color rgb <0.0, 1.0, 1.0>}\n");
            fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
            fprintf(fp,"no_shadow\n");
            fprintf(fp,"}\n");
            fprintf(fp,"\n");

          }
      }
      fprintf(fp,"\n");
      // fprintf(fp,"pigment{color rgb <0.0, 1.0, 1.0>}\n");
    //  fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
    //  fprintf(fp,"no_shadow\n");
      fprintf(fp,"\n");
      //fprintf(fp,"}\n");
      fprintf(fp,"\n");
}
///////////////////////////////////

//------Spheroid visual demo------//
/*object{
   Spheroid( //CenterVector,
             <-1.50,3.00,-2.00>,
             // RadiusVector Rx,Ry,Rz )
             <2.0,1.2,2.5> )
   texture{ pigment{color rgbt<.75,.2,0,.7>}
            finish { phong 1}
          } // end of texture
   scale<1,1,1>
   rotate<0,0,0>
   translate<0,0.0,0>
} //----------------------------------------
*/
///////////////////////////////////
double Spheroid_Cell_povray_visual(){
      fprintf(fp,"object{ Spheroid( <%f,%f,%f>,<%f,%f,%f>)",originx,originy,originz,A_x,A_y,A_z);
      fprintf(fp,"texture{ pigment{ color rgbt< 0.9, 0.9, 0.9, 0.60 > }");
      fprintf(fp,"finish { phong 0.5  phong_size 50}\n");
      fprintf(fp,"}\n");
      fprintf(fp,"finish { reflection 0.0 }\n"); // Glass reflects a bit
      fprintf(fp,"}\n"); //--------------------------------------------------------
}

///////////////////////////////////

/*double Cell_povray_visual(){
      fprintf(fp,"sphere{ <%f,%f,%f>,%f\n",0.0,0.0,0.0,1.0*r_cell);
      fprintf(fp,"texture{ pigment{ color rgbt< 0.9, 0.9, 0.9, 0.60 > }");
      fprintf(fp,"finish { phong 0.5  phong_size 50}\n");
      fprintf(fp,"}\n");
      fprintf(fp,"finish { reflection 0.0 }\n"); // Glass reflects a bit
      fprintf(fp,"}\n"); //--------------------------------------------------------
}*/
//////////////////////////////////
////////////////////////////////////////
double cMT_povray_visual(){

    int i,j;
    //fprintf(fp,"union\n");
    //fprintf(fp,"{\n");
    for(i=0;i<N_cent;i++){
       for(j=0;j<Num_cMT;j++){
          /*if(label_cMT[i][j]==1.0 && l_cmt[i][j]>0.0){
            fprintf(fp,"cylinder{<%f,%f,%f>, <%f,%f,%f>, %f\n",cent_x[i],cent_y[i],cent_z[i],cMTtipx[i][j],cMTtipy[i][j],cMTtipz[i][j],R_cMT_visual);
            fprintf(fp,"\n");
            fprintf(fp,"pigment{color rgb <0.0, 1.0, 0.0>}\n");
            fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
            fprintf(fp,"no_shadow\n");
            fprintf(fp,"\n}");
            //fprintf(fp,"}\n");
            fprintf(fp,"\n");
          }*/
          if(label_cMT[i][j]==-1.0 && l_cmt[i][j]>0.0){
            fprintf(fp,"cylinder{<%f,%f,%f>, <%f,%f,%f>, %f\n",cent_x[i],cent_y[i],cent_z[i],cMTtipx[i][j],cMTtipy[i][j],cMTtipz[i][j],R_cMT_visual);
            fprintf(fp,"\n");
            fprintf(fp,"pigment{color rgb <0.0, 0.0, 1.0>}\n");
            fprintf(fp,"finish{ambient 0.5 specular 0.5 }\n");
            fprintf(fp,"no_shadow\n");
            fprintf(fp,"\n}");
            //fprintf(fp,"}\n");
            fprintf(fp,"\n");
          }
       }
    }

}
////////////////////////////////////////


////////////////////////////////////////////////////
double fancy_visual(int ens,int tstep){
      int i,j,k;

      sprintf(fName,"Ensemble_%d_data_%d.pov",ens,tstep);
      fp = fopen(fName , "w");

      init_povray_visual();

      Chromosome_arm_povray_visual();
      Kinetochore_povray_visual();
      Spheroid_Cell_povray_visual();
      cMT_povray_visual();

      for(i=0;i<N_cent;i++){
         Centrosome_povray_visual(cent_x[i],cent_y[i],cent_z[i]);
      }

      fclose(fp);
}
//////////////////////////////////////////////////////

//.........................................main..................................//

int main(){
    N_ens = 2;
    N_start =1;
    idum = -1255;
    idum2 = -1287;
    A_x = 10.0;
    A_y = 10.0;
    A_z = 10.0;

    for(ens_counter = N_start; ens_counter<N_ens; ens_counter++){
      place_centrosome();
      place_chromosomes();
      Initialize_Arrays();

      t=0;
      while(t<=t_tot){
        Grow_Shrink_cMT();
        rotate_MT();
        cMT_tip();
        Check_cMT_boundary_encounter();
        check_MT_KT_encounter();
        correct_attachments();
        forces();
        move_chromosomes();
        if(t%100 == 0){
          print_chr_kt_coordinates(ens_counter,t);
          fancy_visual(ens_counter,t);
          captured_ch_print(ens_counter,t);
          capture_mt_count(ens_counter,t);
        }
        t++;
      }
    }
}
