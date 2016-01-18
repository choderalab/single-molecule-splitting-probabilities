/**
 * Time-series analysis helper functions for Matlab.
 *
 * @author John D. Chodera.
 */

import java.lang.Exception;

public class timeseries {

  /**
   * Compute transition counts for discrete trajectory.
   *
   * @param A_t   observable trajectory
   * @param max_tau   maximum lag time for accumulating statistics
   * @return the time-correlation function where C_t[tau] is the estimate for lag time tau
   */
  public static double [] compute_correlation_function(double [] A_t, int max_tau) {
    // Determine timeseries length.
    int T = A_t.length;

    // Allocate storage for time-correlation function.
    double [] C_t = new double[max_tau];
    for (int tau = 1; tau <= max_tau; tau++) 
      C_t[tau-1] = 0.0;
    
    // Compute time-correlation function.
    for (int tau = 1; tau <= max_tau; tau++) {
      for(int t0 = 0; t0 < T-tau; t0++)
        C_t[tau-1] += A_t[t0] * A_t[t0+tau];
      C_t[tau-1] /= (double)(T-tau);
    }
    
    return C_t;
  }

  /**
   */
  public static double [] diffusion(double [] x_t, double xbin, double bin_width, int max_tau) {
    // Determine timeseries length.
    int T = x_t.length;

    // Allocate storage for time-correlation function.
    double [] C_t = new double[max_tau];
    for (int tau = 1; tau <= max_tau; tau++) 
      C_t[tau-1] = 0.0;
    
    // Denominator.
    double [] D_t = new double[max_tau];
    for (int tau = 1; tau <= max_tau; tau++) 
      D_t[tau-1] = 0.0;

    // Compute time-correlation function.
    for(int t0 = 0; t0 < T; t0++) {
      if ((xbin-bin_width/2.0 <= x_t[t0]) && (x_t[t0] < xbin+bin_width/2.0)) {
        for (int tau = 1; (tau <= max_tau) && (t0+tau < T); tau++) {
          C_t[tau-1] += (x_t[t0+tau]-xbin)*(x_t[t0+tau]-xbin);
          D_t[tau-1] += 1.0;
        }
      }
    }

    for (int tau = 1; tau <= max_tau; tau++) 
      C_t[tau-1] /= D_t[tau-1];
    
    return C_t;
  }

  /**
   */
  public static double [] diffusion_pande(double [] x_t, double xbin, double bin_width, int max_tau, int tskip) {
    // Determine timeseries length.
    int T = x_t.length;

    // Compute offset.
    double mean = 0.0;
    for (int t = 0; t < T; t++)
      mean += x_t[t];
    double offset = mean;
    //double offset = x_t[0];
    //    for (int t = 0; t < T; t++)
    //      if (x_t[t] < offset)
    //        offset = x_t[t];
    //double offset = 0.0;

    // Allocate storage for time-correlation function.
    double [] C_t = new double[max_tau];
    for (int tau = 1; tau <= max_tau; tau++) 
      C_t[tau-1] = 0.0;
    
    // Denominator.
    double [] D_t = new double[max_tau];
    for (int tau = 1; tau <= max_tau; tau++) 
      D_t[tau-1] = 0.0;

    // Compute time-correlation function.
    for(int t0 = 0; t0 < T; t0 += tskip) {
      if ((xbin-bin_width/2.0 <= x_t[t0]) && (x_t[t0] < xbin+bin_width/2.0)) {
        for (int tau = 1; (tau <= max_tau) && (t0+tau < T); tau++) {
          C_t[tau-1] += (x_t[t0]-offset)*(x_t[t0+tau]-offset);
          D_t[tau-1] += 1.0;
        }
      }
    }
    
    for (int tau = 1; tau <= max_tau; tau++) 
      C_t[tau-1] /= D_t[tau-1];

    return C_t;
  }

  /**
   */
  public static double [] observed_splitting(double [] x_t, double [] bin_edges, double x_A, double x_B) throws Exception {
    // Determine timeseries length.
    int T = x_t.length;

    // Truncate trajectory length so all events end in commitment.
    while ( (x_A < x_t[T-1]) && (x_t[T-1] < x_B) && (T > 0) )
      T--;
    if (T == 0) 
      throw new Exception("No commitment events.");
    
    // Determine number of bins.
    int nbins = bin_edges.length - 1;

    // Accumulate commitment event statistics.
    double [] NA_i = new double[nbins];
    double [] NB_i = new double[nbins];
    int t = 0; // current marker
    int tcommit = 0;
    while (t < T) {
      if ( (x_t[t] <= x_A) || (x_B <= x_t[t]) ) {
        // Reposition current pointer.
        t++;
        tcommit = t;
      } else if ( (x_A < x_t[tcommit]) && (x_t[tcommit] < x_B) ) {
        // Reposition commit pointer.
        tcommit++;
      } else {
        // Determine current bin.
        int i = 0;
        while ((i < nbins) && ((x_t[t] < bin_edges[i]) || (x_t[t] >= bin_edges[i+1])))
          i++;
        // don't accumulate statistics if not within any bin
        if (i >= nbins)
          continue;

        // Determine commitment direction.        
        if (x_t[tcommit] < x_A)
          NA_i[i] += 1.0;
        else
          NB_i[i] += 1.0;                
        
        // Advance current marker.
        t++;
      }
    } 

    // Compute committor to A.
    double [] pA_i = new double[nbins];
    for (int i = 0; i < nbins; i++) {
      if (bin_edges[i+1] <= x_A)
        pA_i[i] = 1.0;
      else if (bin_edges[i] >= x_B)
        pA_i[i] = 0.0;
      else if (NA_i[i] + NB_i[i] > 0.0)
        pA_i[i] = NA_i[i] / (NA_i[i] + NB_i[i]);
    }

    return pA_i;
  }

  /**
   */
  public static double [][] observed_splitting_with_error(double [] x_t, double [] bin_edges, double x_A, double x_B) throws Exception {
    // Determine timeseries length.
    int T = x_t.length;

    // Truncate trajectory length so all events end in commitment.
    while ( (x_A < x_t[T-1]) && (x_t[T-1] < x_B) && (T > 0) )
      T--;
    if (T == 0) 
      throw new Exception("No commitment events.");

    // Determine number of bins.
    int nbins = bin_edges.length - 1;

    // Accumulate commitment event statistics.
    double [] pA_i = new double[nbins];
    double [] dpA_i = new double[nbins];

    for (int i = 0; i < nbins; i++) {
      System.out.printf("%5d / %5d\n", i, nbins);
      if (bin_edges[i+1] <= x_A) {
        pA_i[i] = 1.0;
        dpA_i[i] = 0.0;
      } else if (bin_edges[i] >= x_B) {
        pA_i[i] = 0.0;
        dpA_i[i] = 0.0;
      } else {
        double [] A_t = new double[T];
        double [] B_t = new double[T];

        int t = 0; // current marker
        int tcommit = 0; // commitment pointer
        while (t < T) {
          A_t[t] = 0.0;
          B_t[t] = 0.0;

          if ( (x_t[t] <= x_A) || (x_B <= x_t[t]) ) {  // t outside of [x_A, x_B]
            // Reposition current pointer.
            t++;
            tcommit = t;
          } else if ( (x_A < x_t[tcommit]) && (x_t[tcommit] < x_B) ) { // tcommit inside [x_A, x_B]
            // Reposition commit pointer.
            tcommit++;
          } else {
            // Determine current bin.
            int bin = 0;
            while ((bin < nbins) && ((x_t[t] < bin_edges[bin]) || (x_t[t] >= bin_edges[bin+1])))
              bin++;
            
            if (bin == i) {
              // Determine commitment direction.        
              if (x_t[tcommit] < x_A) {          
                A_t[t] = 1.0;
              } else {
                A_t[t] = 0.0;
              }
              B_t[t] = 1.0;
            }
            
            // Advance current marker.
            t++;
          }
        }

        // Compute correlation times.
        double g_AA = statistical_inefficiency(A_t, A_t);
        double g_AB = statistical_inefficiency(A_t, B_t);
        double g_BB = statistical_inefficiency(B_t, B_t);
        System.out.printf("g_AA = %.1f ; g_AB = %.1f ; g_BB = %.1f\n", g_AA, g_AB, g_BB);

        // Compute means.
        double EA = 0.0;
        double EB = 0.0;
        for (int s = 0; s < T; s++) {
          EA += A_t[s];
          EB += B_t[s];
        }
        EA /= (double)T;
        EB /= (double)T;

        // Compute variances.
        double varA = 0.0;
        double varB = 0.0;
        double covAB = 0.0;
        for (int s = 0; s < T; s++) {
          varA += (A_t[s] - EA)*(A_t[s] - EA);
          varB += (B_t[s] - EB)*(B_t[s] - EB);
          covAB += (A_t[s] - EA)*(B_t[s] - EB);
        }
        varA /= (double)T;
        varB /= (double)T;        
        covAB /= (double)T;        

        // Compute errors.
        double err2A = varA / ((double)T / g_AA);
        double err2B = varB / ((double)T / g_BB);
        double errAB = covAB / ((double)T / g_AB);
        
        // Compute total error.
        pA_i[i] = EA / EB;
        dpA_i[i] = 0.0;
        double err2 = (EA/EB)*(EA/EB)*(err2A/(EA*EA) + err2B/(EB*EB) - 2.0*errAB/(EA*EB));
        if (err2 > 0.0)
          dpA_i[i] = Math.sqrt(err2);
      }
    }

    // Pack return values.
    double [][] retval = new double[2][nbins];
    for (int i = 0; i < nbins; i++) {
      retval[0][i] = pA_i[i];
      retval[1][i] = dpA_i[i];
    }
    
    return retval;
  }

  public static double statistical_inefficiency(double [] A_t, double [] B_t) {
    // Determine timeseries length.
    int T = A_t.length;

    double g = 0.5; // statistical inefficiency

    /* Compute expectation of A and A^2 over trajectory. */
    double EA = 0.0; // Expectation of A over trajectory. 
    double EB = 0.0; // Expectation of B over trajectory. 
    double EAB = 0.0; // Expectation of AB over trajectory. 
    for(int t0 = 0; t0 < T; t0++) {
      EA += A_t[t0];
      EB += B_t[t0];
      EAB += A_t[t0] * B_t[t0];
    }
    EA /= (double)T;
    EB /= (double)T;
    EAB /= (double)T;

    // Return if there is a problem.
    if(EAB - EA*EB == 0)
      return -1.0;

    // Compute unnormalized time-correlation function.
    int tau_step = 1;
    for(int tau = 1; tau < T; tau += tau_step) {      
      //  Compute time-correlation function.
      double C = 0.0;
      
      for(int t0 = 0; t0 < T - tau; t0++) 
	C += (A_t[t0] - EA) * (B_t[t0+tau] - EB) + (B_t[t0] - EB) * (A_t[t0+tau] - EA);

      C /= 2.0 * (double)(T - tau);
      C /= EAB - EA*EB;

      if(C > 0) {
	if(tau == 1)
	  g += C * (1. - (double)tau/T);
	else if(tau == 2)
	  g += C * (1. - (double)tau/T) * 1.5;
	else
	  g += C * (1. - (double)tau/T) * (1 + 0.5*(tau_step-1) + 0.5 * (tau_step));
      }
      else
	break;

      if(tau > 1)
	tau_step++;
    }
  
  g *= 2.0;

  return g;
  }
  
}
