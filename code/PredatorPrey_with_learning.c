

// Within the main() function, various parameters can be modified.
// The time in the code and data files is based on update time, using the time interval as a reference. Thus, 1 update time corresponds to 0.1 time on the graph.

// The random constant generation uses "SFMT.h". For this to work, the "header" and "SFMT.c" files must be present in the same folder as the code.


// Identifier for each individual || 0: DCS chaser, 1: GCS chaser, 2: Target

// Set type representation || '(Number of GCS chasers)*(Total number of chasers + 1) + (Number of DCS chasers)'
// In our simulation, the total number of chasers = 100, so it's equivalent to converting from base 101 to decimal.
// To accommodate scenarios with learning, the maximum set type is considered as '(Total number of chasers + 1)*Total number of chasers'. (assuming all chasers adopt the same strategy).

// The SetNumber file is commented out due to its large size, approximately 200 megabytes. Please refer to it if needed.

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "SFMT.h"

int M, Tmax;   // M : Number of repetitions, Tmax : Total update time during a simulation
int R, NC, NCd, NCg, NT;  // R : Disk radius, N : Number of individuals
int ini_NC, ini_NCd, ini_NCg, ini_NT;  // Number of initial individuals
int r_danger, r_learn;   // r_danger(=r_haz) : Maximal distance for recognizing hazard, r_learn : Maximal distance for learning
int ht, rt;   // ht : Maximum amount of update for the hunting mode, hr : Maximum amount of update for the rest mode
int catch_count;   // catch_count : Number of hunts (for data file recording)
double vdcdt, vgcdt, vtdt;  // Distance moved per single update
double r_space;  // (=r_min) Minimal distance between chaser or targets
int dt;   // Time interval for writing data file
int CDCS, CGCS, SDCS, SGCS;   // CDCS : Number of targets caught by DCS in the hetero set, CGCS : Number of targets caught by GCS in the hetero set, SDCS : Number of targets caught by DCS in the homo set, SGCS : Number of targets caught by GCS in the homo set
int beforefix;  // Check if fixation has occurred

double **C, **T;   // C : Position of chasers, T : Position of targets
double **Ctemp, **Ttemp;   // Ctemp : Temporary position of chasers (for the update process), Ttemp : Temporary position of targets (for the update process)
int *CS, *CS_temp;  // CS : Strategy of chasers,  CS_temp : for learning process
int *CC, *CP;   // CC : Status of the chaser (hunting or rest), CP : Probability of the chaser's status
double **VC, **VT;   // VC : Velocity of chasers, VT : Velocity of targets
double **Dis, *LearnDis;  // Dis : Distance between chasers and targets, LearnDis : Distance between chasers when learning occurs
int *Set, *Set_temp, **CatchSet;   // Set: The set chasing the target, Set_temp: For cases where multiple targets are caught in a single update, CatchSet: The set that successfully hunted the target
int *C_order, *T_order, *NET, *NEC;  // order : Movement order, NET : Nearest target, NEC : Nearest chaser
int **CatchNumber, **BelongSet, **SetNumber, **Number, **FinishNumber, **Fixation;   // Arrays for data storage; for more details, refer to the "write_file()" function.


int m,t;
int seed;
sfmt_t sfmt;

void main_setting(), ini_setting(), main_free(), ini_free();
void calculate_distance(), find_nearest(), check_set(), chaser_velocity(), target_velocity(), chaser_move(), target_move(), record_catch_number(), check_fixation(), record_number(), no_fixation(), record_finish_number();
void write_file_Position(), write_file_CatchSet(), write_file_CatchNumber(), write_file_BelongSet(),write_file_SetNumber(), write_file_Number(), write_file_FinishNumber(), write_file_Fixation();



int main(){

  M = 2; Tmax = 1000000;
  ini_NCd = 50; ini_NCg = 50; ini_NT = 50;
  vdcdt = 0.082; vgcdt = 0.070; vtdt = 0.10;
  ini_NC = ini_NCd+ini_NCg;
  
  R = 500; 
  r_danger = 50, r_learn = 200;
  ht = 1000; rt = 1000;   
  catch_count = 0;
  r_space = 0.15;
  dt = 1000;
  
  main_setting();
  for(m=0;m<M;m++){

    seed = time(NULL);
    sfmt_init_gen_rand(&sfmt,seed);
    
    NC = ini_NC; NCd = ini_NCd; NCg = ini_NCg; NT = ini_NT;
    CDCS = 0; CGCS = 0; SDCS = 0; SGCS = 0;
    t = 0;
    beforefix = 1;
    ini_setting();
    for(t=0;t<Tmax;t++){
      if(m==0 && t%100000==0) write_file_Position();
      calculate_distance();
      find_nearest();
      check_set();
      chaser_velocity();
      target_velocity();
      chaser_move();
      target_move();
      if(beforefix) check_fixation();
      if(t%dt==0){
	record_catch_number();
	record_number();
      }
    }
    if(beforefix) no_fixation();
    record_finish_number();
    ini_free();
  }
  write_file_CatchSet();
  write_file_CatchNumber();
  write_file_BelongSet();
  // write_file_SetNumber();
  write_file_Number();
  write_file_FinishNumber();
  write_file_Fixation();
  main_free();
} 











// ********************** Base setting and calculation code ***********************


// a: me, b: opponent, return: where the opponent is in relation to me                                           
double delta(double a, double b){
  return b-a;
}


double distance(double x, double y){
  return sqrt(x*x+y*y);
}


void ini_setting(){
  C = (double**)calloc(ini_NC,sizeof(double*));
  T = (double**)calloc(ini_NT,sizeof(double*));
  Ctemp = (double**)calloc(ini_NC,sizeof(double*));
  Ttemp = (double**)calloc(ini_NT,sizeof(double*));
  VC = (double**)calloc(ini_NC,sizeof(double*));
  VT = (double**)calloc(ini_NT,sizeof(double*));
  Dis = (double**)calloc(ini_NC,sizeof(double*));
  for(int i=0;i<ini_NC;i++){
    C[i] = (double*)calloc(2,sizeof(double));
    Ctemp[i] = (double*)calloc(3,sizeof(double));
    VC[i] = (double*)calloc(2,sizeof(double));
    Dis[i] = (double*)calloc(ini_NT,sizeof(double));
  }
  // chaser's initial position (consider r_space)                                           
  for(int i=0;i<ini_NC;i++){
    double r = sfmt_genrand_res53(&sfmt)*R;
    double theta = sfmt_genrand_res53(&sfmt)*2.*M_PI;
    C[i][0] = r*cos(theta);
    C[i][1] = r*sin(theta);
    int n = 0;
    for(int j=0;j<i;j++){
      if(distance(delta(C[j][0],C[i][0]),delta(C[j][1],C[i][1])) < r_space) break;
      n += 1;
    }
    if(n!=i) i -= 1;
  }
    
  for(int i=0;i<ini_NT;i++){
    T[i] = (double*)calloc(2,sizeof(double));
    Ttemp[i] = (double*)calloc(3,sizeof(double));
    VT[i] = (double*)calloc(2,sizeof(double));
  }
  // target's initial position (consider r_space)                                           
  for(int i=0;i<ini_NT;i++){
    double r = sfmt_genrand_res53(&sfmt)*R;
    double theta = sfmt_genrand_res53(&sfmt)*2.*M_PI;
    T[i][0] = r*cos(theta);
    T[i][1] = r*sin(theta);
    int n = 0;
    for(int j=0;j<i;j++){
      if(distance(delta(T[j][0],T[i][0]),delta(T[j][1],T[i][1])) < r_space) break;
      n += 1;
    }
    if(n!=i) i -= 1;
  }
  
  CS = (int*)calloc(ini_NC,sizeof(int));
  CS_temp = (int*)calloc(ini_NC,sizeof(int));
  CC = (int*)calloc(ini_NC,sizeof(int));
  CP = (int*)calloc(ini_NC,sizeof(int));
  for(int i=ini_NCd;i<ini_NC;i++){
    CS[i] = 1;
    CS_temp[i] = 1;
  }
  for(int i=0;i<ini_NC;i++) CC[i] = 1;
  
  Set = (int*)calloc(ini_NT,sizeof(int));
  Set_temp = (int*)calloc(ini_NT,sizeof(int));
}


void main_setting(){
  C_order = (int*)calloc(ini_NC,sizeof(int));
  T_order = (int*)calloc(ini_NT,sizeof(int));
  NET = (int*)calloc(ini_NC,sizeof(int));
  NEC = (int*)calloc(ini_NT,sizeof(int));
  LearnDis = (double*)calloc(M*ini_NT*1000,sizeof(double));
  
  CatchSet = (int**)calloc(M*ini_NT*1000,sizeof(int*));
  for(int i=0;i<M*ini_NT*1000;i++) CatchSet[i] = (int*)calloc(4+ini_NT,sizeof(int));

  CatchNumber = (int**)calloc(Tmax/dt,sizeof(int*));
  for(int i=0;i<Tmax/dt;i++) CatchNumber[i] = (int*)calloc(5,sizeof(int));

  BelongSet = (int**)calloc(Tmax/dt,sizeof(int*));
  for(int i=0;i<Tmax/dt;i++) BelongSet[i] = (int*)calloc(4,sizeof(int));

  SetNumber = (int**)calloc(Tmax/dt,sizeof(int*));
  for(int i=0;i<Tmax/dt;i++) SetNumber[i] = (int*)calloc((ini_NC+1)*(ini_NC+1),sizeof(int));

  Number = (int**)calloc(Tmax/dt,sizeof(int*));
  for(int i=0;i<Tmax/dt;i++) Number[i] = (int*)calloc(2,sizeof(int));

  FinishNumber = (int**)calloc(M,sizeof(int*));
  for(int i=0;i<M;i++) FinishNumber[i] = (int*)calloc(3,sizeof(int));

  Fixation = (int**)calloc(M,sizeof(int*));
  for(int i=0;i<M;i++) Fixation[i] = (int*)calloc(3,sizeof(int));
  
}


void ini_free(){
  for(int i=0;i<ini_NC;i++){
    free(C[i]);
    free(Ctemp[i]);
    free(VC[i]);
    free(Dis[i]);
  } 
  for(int i=0;i<ini_NT;i++){
    free(T[i]);
    free(Ttemp[i]);
    free(VT[i]);
  }
  free(C);
  free(Ctemp);
  free(VC);
  free(Dis);
  free(T);
  free(Ttemp);
  free(VT);
  free(CS);
  free(CS_temp);
  free(CC);
  free(CP);
  free(Set);
  free(Set_temp);
}


void main_free(){
  free(C_order);
  free(T_order);
  free(NET);
  free(NEC);
  
  for(int i=0;i<M*ini_NT*1000;i++) free(CatchSet[i]);
  free(CatchSet);

  for(int i=0;i<Tmax/dt;i++) free(CatchNumber[i]);
  free(CatchNumber);

  for(int i=0;i<Tmax/dt;i++) free(BelongSet[i]);
  free(BelongSet);

  for(int i=0;i<Tmax/dt;i++) free(SetNumber[i]);
  free(SetNumber);

  for(int i=0;i<Tmax/dt;i++) free(Number[i]);
  free(Number);

  for(int i=0;i<M;i++){
    free(FinishNumber[i]);
    free(Fixation[i]);
  }
  free(FinishNumber);
  free(Fixation);
}


// If crossing the (x,y) boundary, move to the edge.
void check_boundary(double *x,double *y){   
  double d = distance(*x,*y);
  if(d>R){
    double theta = asin(*y/d);
    if(*x<0){
      *x = -R*cos(theta);
      *y = R*sin(theta);
    }
    else{
      *x = R*cos(theta);
      *y = R*sin(theta);
    }
  }
}


// Return the distance from (x,y) to the closest boundary.
double boundary_distance(double x, double y){
  return R - distance(x,y);
}


// Generate a random order for the chaser.
void C_order_setting(){
  int a,temp;
  for(int i=0;i<NC;i++) C_order[i] = i;
  for(int i=0;i<NC;i++){
    a = rand()%NC;
    temp = C_order[a];
    C_order[a] = C_order[i];
    C_order[i] = temp;
  }
}


// Generate a random order for the target.
void T_order_setting(){
  int a,temp;
  for(int i=0;i<NT;i++) T_order[i] = i;
  for(int i=0;i<NT;i++){
    a =	rand()%NT;
    temp = T_order[a];
    T_order[a] = T_order[i];
    T_order[i] = temp;
  }
}



// ********************** Simulation code *********************** 


// Calculating the distance between the chaser and the target -> storing in an array
void calculate_distance(){
  for(int ci=0;ci<NC;ci++){
    for(int ti=0;ti<NT;ti++){
      Dis[ci][ti] = distance(delta(C[ci][0],T[ti][0]),delta(C[ci][1],T[ti][1]));
    }
  }
}


// Finding the nearest for each individual
void find_nearest(){
  for(int ci=0;ci<NC;ci++) NET[ci] = 0;
  for(int ti=0;ti<NT;ti++){
    NEC[ti] = 0;
    for(int ci=0;ci<NC;ci++){
      if(Dis[ci][ti] < Dis[ci][NET[ci]]) NET[ci] = ti;   // The closest ti to ci
      if(Dis[ci][ti] < Dis[NEC[ti]][ti]) NEC[ti] = ci;   // The closest ci to ti
    }
    if(boundary_distance(T[ti][0],T[ti][1]) < Dis[NEC[ti]][ti]) NEC[ti] = 1000;  // If the boundary is closest to the target, set it to 1000
  }
}


// Check and store the composition of the set in an array
void check_set(){
  for(int ti=0;ti<ini_NT;ti++) Set[ti] = 0;
  for(int ci=0;ci<NC;ci++){
    int ti = NET[ci];
    if(CS[ci]) Set[ti] += 1;
    else Set[ti] += ini_NC+1;
  }
  for(int ti=0;ti<ini_NT;ti++){
    Set_temp[ti] = Set[ti];
    SetNumber[t/dt][Set[ti]]++;
  }
    
  for(int ti=0;ti<NT;ti++){
    int d = Set[ti]/(ini_NC+1);
    int g = Set[ti]%(ini_NC+1);
    if(d!=0 && g!=0){
      BelongSet[t/dt][0] += d;
      BelongSet[t/dt][1] += g;
    }
    else if(g==0) BelongSet[t/dt][2] += d;
    else BelongSet[t/dt][3] += g;
  }
}


// Velocity of chasers
void chaser_velocity(){
  for(int ci=0;ci<NC;ci++){
    
    double v;
    if(CS[ci]) v = vgcdt;
    else v = vdcdt;
    
    // hunting
    if(CC[ci]){
      int ti = NET[ci];

      // catch
      if(Dis[ci][ti]<v){
	VC[ci][0] = delta(C[ci][0],T[ti][0]);
	VC[ci][1] = delta(C[ci][1],T[ti][1]);
      }

      // GCS
      else if(CS[ci]){
	double xsum = 0;
	double ysum = 0;
	int n = 0;
	for(int cj=0;cj<NC;cj++){
	  if(cj!=ci && Dis[cj][ti]<=Dis[ci][ti]){
	    xsum += C[cj][0];
	    ysum += C[cj][1];
	    n += 1;
	  }
	}
	double cx = (n+1)*T[ti][0] - xsum;
	double cy = (n+1)*T[ti][1] - ysum;
	double dx = delta(C[ci][0],cx);
	double dy = delta(C[ci][1],cy);
	double d = distance(dx,dy);
	VC[ci][0] = v*dx/d;
	VC[ci][1] = v*dy/d;
      }

      // DCS
      else{
	double dx = delta(C[ci][0],T[ti][0]);
	double dy = delta(C[ci][1],T[ti][1]);
	double d = distance(dx,dy);
	VC[ci][0] = v*dx/d;
	VC[ci][1] = v*dy/d;
      }
    }
    
    // rest
    else{
      double theta = sfmt_genrand_res53(&sfmt)*6.28;
      VC[ci][0] = v*cos(theta);
      VC[ci][1] = v*sin(theta);
    }
  }
}
 

// Target's death and reborn
void target_death(int ti,int ci){
  CatchSet[catch_count][0] = m;
  CatchSet[catch_count][1] = t;
  CatchSet[catch_count][2] = Set[ti];  
  CatchSet[catch_count][3] = CS[ci];
  for(int tj=0;tj<ini_NT;tj++){
    CatchSet[catch_count][4+tj] = Set_temp[tj];
  }    
  catch_count++;

  if(Set[ti]<=ini_NC) SGCS++;    // GCS in homo set
  else if(Set[ti]%(ini_NC+1)==0) SDCS++;   // DCS in homo set
  else if(CS[ci]) CGCS++;   // GCS in hetero set
  else CDCS++;   // DCS in hetero set

  for(int cj=0;cj<NC;cj++){
    if(Dis[cj][ti] <= r_learn) CS_temp[cj] = CS_temp[ci];
  }

  // target's reborn
  int space = 1;
  while(space){
    double r = sfmt_genrand_res53(&sfmt)*R;
    double theta = sfmt_genrand_res53(&sfmt)*2.*M_PI;
    T[ti][0] = r*cos(theta);
    T[ti][1] = r*sin(theta);
    int n = 0;
    // for r_space                                                                          
	 for(int i=0;i<ini_NT;i++){
	   if(distance(delta(T[i][0],T[ti][0]),delta(T[i][1],T[ti][1])) > r_space) n += 1;
	 }
    if(n==ini_NT-1) space = 0;
  }
}


// Change in the chaser's state (hunting or rest)
void change_condition(){
  for(int ci=0;ci<NC;ci++){
    CP[ci]++;
    double p = sfmt_genrand_res53(&sfmt);
    if(CC[ci] && (double)CP[ci]/ht >= p){  // hunt -> rest
      CC[ci] = 0;  
      CP[ci] = 0;
    }
    else if((double)CP[ci]/ht >= p){  // rest -> hunt
      CC[ci] = 1;   
      CP[ci] = 0;
    }
  }
}


// movement of chasers
void chaser_move(){

  change_condition();
  
  C_order_setting();
  for(int i=0;i<NC;i++){
    int ci = C_order[i];
    int ti = NET[ci];

    Ctemp[i][0] = C[ci][0] + VC[ci][0];
    Ctemp[i][1] = C[ci][1] + VC[ci][1];
    Ctemp[i][2] = ci;
    check_boundary(&Ctemp[i][0],&Ctemp[i][1]);
    
    // check r_space
    int n = 0;
    for(int cj=0;cj<i;cj++){
      if(distance(delta(Ctemp[i][0],Ctemp[cj][0]),delta(Ctemp[i][1],Ctemp[cj][1])) < r_space) break;
      n += 1;
    }
    // If space condition is not met, return to the original position
    if(n!=i){
      Ctemp[i][0] = C[ci][0];
      Ctemp[i][1] = C[ci][1];
    }
    if(Ctemp[i][0]==T[ti][0] && Ctemp[i][1]==T[ti][1]) target_death(ti,ci);
  }

  // Save the new position in an array
  for(int i=0;i<NC;i++){
    int ci = (int)Ctemp[i][2];
    C[ci][0] = Ctemp[i][0];
    C[ci][1] = Ctemp[i][1];
  }
  NCd = 0;
  NCg = 0;
  for(int ci=0;ci<NC;ci++){
    CS[ci] = CS_temp[ci];   // change the strategy by learning
    if(CS[ci]) NCg++;
    else NCd++;
  }
}  


// velocoty of targets
void target_velocity(){
  for(int ti=0;ti<NT;ti++){
    int ci = NEC[ti];
    // boundary
    if(NEC[ti]==1000){
      // escape
      if(boundary_distance(T[ti][0],T[ti][1]) < r_danger){
	double dx = -T[ti][0];
	double dy = -T[ti][1];
	double d = distance(dx,dy);
	VT[ti][0] = vtdt*dx/d;
	VT[ti][1] = vtdt*dy/d;
      }
      // random
      else{
	double theta = sfmt_genrand_res53(&sfmt)*6.28;
	VT[ti][0] = vtdt*cos(theta);
        VT[ti][1] = vtdt*sin(theta);
      }	
    }
    // chaser
    else{
      // escape
      if(Dis[ci][ti] < r_danger){
	double dx = delta(C[ci][0],T[ti][0]);
	double dy = delta(C[ci][1],T[ti][1]);
	double d = Dis[ci][ti];
	VT[ti][0] = vtdt*dx/d;
	VT[ti][1] = vtdt*dy/d;
      }
      else{
	// random
	double theta = sfmt_genrand_res53(&sfmt)*6.28;
        VT[ti][0] = vtdt*cos(theta);
        VT[ti][1] = vtdt*sin(theta);
      }
    }
  }
}


// Movement of targets
void target_move(){
  
  T_order_setting();
  for(int i=0;i<NT;i++){
    int ti = T_order[i];    
    Ttemp[i][0] = T[ti][0] + VT[ti][0];
    Ttemp[i][1] = T[ti][1] + VT[ti][1];
    Ttemp[i][2] = ti;
    check_boundary(&Ttemp[i][0],&Ttemp[i][1]);

    // check r_space
    int n = 0;
    for(int tj=0;tj<i;tj++){
      if(distance(delta(Ttemp[i][0],Ttemp[tj][0]),delta(Ttemp[i][1],Ttemp[tj][1])) < r_space) break;
      n += 1;
    }
    // If space condition is not met, return to the original position
    if(n!=i){
      Ttemp[i][0] = T[ti][0];
      Ttemp[i][1] = T[ti][1];
    }
  }

  // Save the new position in an array
  for(int i=0;i<NT;i++){
    int ti = (int)Ttemp[i][2];
    T[ti][0] = Ttemp[i][0];
    T[ti][1] = Ttemp[i][1];
  }
}


void record_catch_number(){
  CatchNumber[t/dt][0] += CDCS;
  CatchNumber[t/dt][1] += CGCS;
  CatchNumber[t/dt][2] += SDCS;
  CatchNumber[t/dt][3] += SGCS;
  CatchNumber[t/dt][4] += CDCS+CGCS+SDCS+SGCS;
}


void record_number(){
  Number[t/dt][0] += NCd;
  Number[t/dt][1] += NCg;
}


void record_finish_number(){
  FinishNumber[m][0] = NCd;
  FinishNumber[m][1] = NCg;
  FinishNumber[m][2] = CDCS+CGCS+SDCS+SGCS;
}


void check_fixation(){
  if(NCd==ini_NC){
    Fixation[m][0] = t;
    Fixation[m][2] = CDCS+CGCS+SDCS+SGCS;
    beforefix = 0;
  }
  if(NCg==ini_NC){
    Fixation[m][0] = t;
    Fixation[m][1] = 1;
    Fixation[m][2] = CDCS+CGCS+SDCS+SGCS;
    beforefix = 0;
  }
}


void no_fixation(){
  Fixation[m][0] = t;
  Fixation[m][1] = -1;
  Fixation[m][2] = CDCS+CGCS+SDCS+SGCS;
}

 
// *************************** Data file code ****************************


// Locations of individual according to updates (adjust time interval as desired in the main function, recorded only for a single execution)
// col1: x position, col2: y position, col3: type of individual(0:DCS, 1:GCS, 2:Target) 
void write_file_Position(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/01_Position_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf_t%d.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space,t);
  wf = fopen(fn,"w");
  fprintf(wf,"# Locations of individual according to updates (adjust time interval as desired in the main function, recorded only for a single execution)\n");
  fprintf(wf,"# col1: x position, col2: y position, col3: type of individual(0:DCS, 1:GCS, 2:Target)\n");
  for(int ci=0;ci<NC;ci++) fprintf(wf,"%lf %lf %d\n",C[ci][0],C[ci][1],CS[ci]);
  for(int ti=0;ti<NT;ti++) fprintf(wf,"%lf %lf %d\n",T[ti][0],T[ti][1],2);
  fclose(wf);
}
    

// Set that caught the target
// col1: Simulation number, col2: Update time, col3: Set type, col4: Chaser strategy for catching the target, col5 to col54: Set types chasing each target
void write_file_CatchSet(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/02_Catch_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Set that caught the target\n");
  fprintf(wf,"# col1: Simulation number, col2: Update time, col3: Set type, col4: Chaser strategy for catching the target, col5 to col54: Set types chasing each target\n");
  for(int i=0;i<catch_count;i++){
    for(int j=0;j<4+ini_NT;j++) fprintf(wf,"%d ",CatchSet[i][j]);
    fprintf(wf,"\n");
  }
  fclose(wf);
}


// Number of caught targets in hetero or homo sets (accumulated over time)
// col1: Update time, col2: Number of targets caught by DCS in hetero sets, col3: Number of targets caught by GCS in hetero sets, col4: Number of targets caught by DCS in homo sets, col5: Number of targets caught by GCS in homo sets, col6: Total number of caught targets
void write_file_CatchNumber(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/03_CatchNumber_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Number of caught targets in hetero or homo sets (accumulated over time)\n");
  fprintf(wf,"# col1: Update time, col2: Number of targets caught by DCS in hetero sets, col3: Number of targets caught by GCS in hetero sets, col4: Number of targets caught by DCS in homo sets, col5: Number of targets caught by GCS in homo sets, col6: Total number of caught targets\n");
  for(int tt=0;tt<Tmax/dt;tt++) fprintf(wf,"%d %lf %lf %lf %lf %lf\n",tt*dt,(double)CatchNumber[tt][0]/M,(double)CatchNumber[tt][1]/M,(double)CatchNumber[tt][2]/M,(double)CatchNumber[tt][3]/M,(double)CatchNumber[tt][4]/M);
  fclose(wf);
}
 


// Number of chasers in hetero or homo sets (within a specific time interval)
// col1: Update time, col2: Number of DCS chasers in hetero sets, col3: Number of GCS chasers in hetero sets, col4: Number of DCS chasers in homo sets, col5: Number of GCS chasers in homo sets
void write_file_BelongSet(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/04_BelongSet_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Number of chasers in hetero or homo sets (within a specific time interval)\n");
  fprintf(wf,"# col1: Update time, col2: Number of DCS chasers in hetero sets, col3: Number of GCS chasers in hetero sets, col4: Number of DCS chasers in homo sets, col5: Number of GCS chasers in homo sets\n");
  for(int tt=0;tt<Tmax/dt;tt++) fprintf(wf,"%d %lf %lf %lf %lf\n",tt*dt,(double)BelongSet[tt][0]/M,(double)BelongSet[tt][1]/M,(double)BelongSet[tt][2]/M,(double)BelongSet[tt][3]/M);
  fclose(wf);
}


// Number of sets belonging to the corresponding set type (within a specific time interval)
// col1: Update time, col2 to col10202(for 100 chasers): Number of sets belonging to the respective type
void write_file_SetNumber(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/05_SetNumber_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Number of sets belonging to the corresponding set type (within a specific time interval)\n");
  fprintf(wf,"# col1: Update time, col2 to col10202(for 100 chasers): Number of sets belonging to the respective type\n");
  for(int tt=0;tt<Tmax/dt;tt++){
    fprintf(wf,"%d ",tt*dt);
    for(int j=0;j<(ini_NC+1)*(ini_NC+1);j++) fprintf(wf,"%lf ",(double)SetNumber[tt][j]/M);
    fprintf(wf,"\n");
  }
  fclose(wf);
}


// Number of chasers based on strategy (for the given time)
// col1: Update time, col2: Number of DCS, col3: Number of GCS
void write_file_Number(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/06_Number_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Number of chasers based on strategy (for the given time)\n");
  fprintf(wf,"# col1: Update time, col2: Number of DCS, col3: Number of GCS\n");
  for(int tt=0;tt<Tmax/dt;tt++) fprintf(wf,"%d %lf %lf\n",tt*dt,(double)Number[tt][0]/M,(double)Number[tt][1]/M);
  fclose(wf);
}


// Number of individuals at the end of the simulation
// col1: Simulation number, col2: Number of DCS chasers at the end of the simulation, col3: Number of GCS chasers at the end of the simulation, col4: Number of targets captured during that simulation
void write_file_FinishNumber(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/07_FinishNumber_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Number of individuals at the end of the simulation\n");
  fprintf(wf,"# col1: Simulation number, col2: Number of DCS chasers at the end of the simulation, col3: Number of GCS chasers at the end of the simulation, col4: Number of targets captured during that simulation\n");
  for(int mm=0;mm<M;mm++) fprintf(wf,"%d %d %d %d\n",mm,FinishNumber[mm][0],FinishNumber[mm][1],FinishNumber[mm][2]);
  fclose(wf);
}


// Information on fixation for each simulation
// col1: Simulation number, col2: fixation time, col3: fixation strategy (0:DCS, 1:GCS, -1:no fixation occurred), col4: Number of targets captured up to the point of fixation
void write_file_Fixation(){
  FILE *wf;
  char fn[300];
  sprintf(fn,"../result/08_Fixation_withLearning_M%d_Tmax%d_R%d_NCd%d_NCg%d_NT%d_vdcdt%.3lf_vgcdt%.3lf_vtdt%.3lf_rd%d_rl%d_ht%d_rs%.2lf.txt",M,Tmax,R,ini_NCd,ini_NCg,ini_NT,vdcdt,vgcdt,vtdt,r_danger,r_learn,ht,r_space);
  wf = fopen(fn,"w");
  fprintf(wf,"# Information on fixation for each simulation\n");
  fprintf(wf,"# col1: Simulation number, col2: Fixation time, col3: Fixation strategy (0:DCS, 1:GCS, -1:no fixation occurred), col4: Number of targets captured up to the point of fixation\n");
  for(int mm=0;mm<M;mm++) fprintf(wf,"%d %d %d %d\n",mm,Fixation[mm][0],Fixation[mm][1],Fixation[mm][2]);
  fclose(wf);
}


 














 
