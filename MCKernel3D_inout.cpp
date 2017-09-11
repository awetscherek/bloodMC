
#include "mex.h"

#include <random>
#include <time.h>

// x0, y0, z0   : initial coordinates
// L1, L2       : width = depth and height of unit cell 
// Lr, Lz       : radius and HALF height of cylindical RBC
// dr_e, dr_p   : 2 * step size (inside and outside of RBC)
// k2D_e, k2D_p : probability measure for membrane transit 
//                (coming from inside or outside of RBC)
// NSteps       : number of steps for random walk
void MCKernel3D_inout(double x0, double y0, double z0, 
                      double L1, double L2, double Lr, double Lz, 
                      double dr_e, double dr_p, double k2D_e, double k2D_p, 
                      int NSteps, double *outmatrix, mwSize nn) {
    
    double x, y, z;        // current coordinates
    double r_cell2;        // squared radial coordinate in unit cell
    double z_cell;         // z-coordinate in unit cell
    double ds, ds_1, ds_2; // for distance to RBC membrane
    int transit;           // stores whether membrane transit is happening
                           // 0: no transit, 1: forbidden, 2: transit
    
    // initialize rng - uniformly distributed between -0.5 and +0.5
    std::mt19937 generator((unsigned int) time(0) 
       + (int) (fabs(x0) / L1 * 1E6) + (int) (fabs(y0) / L1 * 1E6));
    std::uniform_real_distribution<double> dis(-0.5, 0.5);
    
    // variables for summing up center of mass
    double xcm = 0.0, ycm = 0.0, zcm = 0.0;
    
    // time spend within RBC, outside RBC, time of first membrane transit
    double t_in = 0.0, t_out = 0.0, t_first = 0.0;
    
    // current x- and y-coordinate relative to unit cell
    double x_cell = x0 - floor(x0 / L1) * L1 - L1 / 2.0;  
    double y_cell = y0 - floor(y0 / L1) * L1 - L1 / 2.0;
    
    // initial z-coordinate relative to unit cell
    double z0_cell  = z0 - floor(z0 / L2) * L2 - L2 / 2.0;
    
    // squared radial coordinate  and z-coordinate relative to unit cell
    double r0_cell2 = x_cell * x_cell + y_cell * y_cell;
    double z_cell2  = z0_cell * z0_cell;                  
    
    double Lr2 = Lr * Lr; // squared radius of RBC cylinder
    double Lz2 = Lz * Lz; // squared HALF height of RBC cylinder
   
     // Determine whether starting location is in RBC: 
    int x0_in_RBC = ((r0_cell2 < Lr2) && (z_cell2 < Lz2));
     
    for (int cStep = 0; cStep < NSteps; cStep++) { // loop over timesteps
        
        do { // repeat as long as step is forbidden ...
           
            if (x0_in_RBC) { // inside: use dr_e and not dr_p ... 
                
                // update coordinates
                x = x0 + dr_e * dis(generator);
                y = y0 + dr_e * dis(generator);
                z = z0 + dr_e * dis(generator);
                x_cell = x - floor(x / L1) * L1 - L1 / 2.0;
                y_cell = y - floor(y / L1) * L1 - L1 / 2.0;
                z_cell = z - floor(z / L2) * L2 - L2 / 2.0;
                r_cell2 = x_cell * x_cell + y_cell * y_cell;
                z_cell2 = z_cell * z_cell;
                
                // check if transit happens
                transit = ((r_cell2 >= Lr2) || (z_cell2 >= Lz2));
                
                if (transit) {
                    
                    // calculate min. distance to cell membrane
                    ds_1 = Lz - z0_cell;    
                    ds_2 = Lz + z0_cell; 
                    ds   = Lr - sqrt(r0_cell2);
                    if (ds_1 < ds) ds = ds_1;
                    if (ds_2 < ds) ds = ds_2;
                    
                    // permeability determines if transit is allowed:
                    transit += (dis(generator) + 0.5 < ds * k2D_e);       
                    
                    // if transit is allowed
                    if (transit == 2) {
                           
                        // not outside before => it is the first transit!
                        if (t_out < 1e-12) t_first = (double) cStep + 0.5; 
                        
                        // update residence times:
                        t_in = (double) cStep + 0.5 - t_out;
                    }
                }
                    
            } else {
                
                // update coordinates
                x = x0 + dr_p * dis(generator);
                y = y0 + dr_p * dis(generator);
                z = z0 + dr_p * dis(generator);
                x_cell = x - floor(x / L1) * L1 - L1 / 2.0;
                y_cell = y - floor(y / L1) * L1 - L1 / 2.0;
                z_cell = z - floor(z / L2) * L2 - L2 / 2.0;
                r_cell2 = x_cell * x_cell + y_cell * y_cell;
                z_cell2 = z_cell * z_cell;
                
                // check if transit happens
                transit = ((r_cell2 < Lr2) && (z_cell2 < Lz2));
                
                if (transit) { 
                    
                    // calculate min. distance to cell membrane
                    ds_1 = fabs(Lz - z0_cell);
                    ds_2 = fabs(Lz + z0_cell);
                    ds   = fabs(Lr - sqrt(r0_cell2));
                    if (ds_1 < ds) ds = ds_1;
                    if (ds_2 < ds) ds = ds_2;
                    
                    // permeability determines if transit is allowed:
                    transit += (dis(generator) + 0.5 < ds * k2D_p);
                    
                    // if transit is allowed
                    if (transit == 2) {
                        
                        // not inside before => it is the first transit!
                        // negative to distinguish from starting inside
                        if (t_in < 1e-12) t_first = -((double) cStep + 0.5); 
                        
                        // update residence times:
                        t_out = (double) cStep + 0.5 - t_in;
                    }
                }
            }
            
        } while (transit == 1);
        
        // add coordinates to center of mass:
        xcm += x;
        ycm += y;
        zcm += z;
        
        // update coordinates:
        x0 = x;
        y0 = y;
        z0 = z;
        r0_cell2 = r_cell2;
        z0_cell = z_cell;
               
        // change status, if transit has happened in this step
        if (transit == 2) x0_in_RBC = 1 - x0_in_RBC;
    }
            
    // export center of mass coordinates:
    outmatrix[0] = xcm / NSteps;
    outmatrix[1] = ycm / NSteps;  
    outmatrix[2] = zcm / NSteps;
        
    // export residence times:
    if (x0_in_RBC > 0.5) {
       outmatrix[3] = (double) NSteps - t_out;
       outmatrix[4] = t_out;
    } else {
       outmatrix[3] = t_in;
       outmatrix[4] = (double) NSteps - t_in;
    }
    
    // export first transit time:
    outmatrix[5] = t_first;
    
    // export final coordinates (for next interval)
    outmatrix[6] = x0;
    outmatrix[7] = y0;
    outmatrix[8] = z0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
  if (nrhs != 12) mexErrMsgIdAndTxt("MCKernel3D_inout", 
    "12 inputs required. x0, y0, z0, L1, L2, Lr, Lz, dr_e, dr_p, k2D_e, k2D_p, NSteps.");
  if (nlhs !=  1) mexErrMsgIdAndTxt("MCKernel3D_inout", 
    "One output required.");
    
  double x0	       = mxGetScalar(prhs[ 0]);
  double y0	       = mxGetScalar(prhs[ 1]);
  double z0	       = mxGetScalar(prhs[ 2]);
  double L1        = mxGetScalar(prhs[ 3]);
  double L2        = mxGetScalar(prhs[ 4]);
  double Lr        = mxGetScalar(prhs[ 5]);
  double Lz        = mxGetScalar(prhs[ 6]);
  double dr_e      = mxGetScalar(prhs[ 7]);
  double dr_p      = mxGetScalar(prhs[ 8]);
  double k2D_e     = mxGetScalar(prhs[ 9]);
  double k2D_p     = mxGetScalar(prhs[10]);
  int NSteps       = static_cast<int>(mxGetScalar(prhs[11]) + 0.5);
       
  mwSize ncols = 9; 
  plhs[0] = mxCreateDoubleMatrix(1, ncols, mxREAL);

  MCKernel3D_inout(x0, y0, z0, L1, L2, Lr, Lz, dr_e, dr_p,
                   k2D_e, k2D_p, NSteps, mxGetPr(plhs[0]), ncols);
}
