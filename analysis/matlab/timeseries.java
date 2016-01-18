/**
 * Time-series analysis helper functions for Matlab.
 *
 * @author John D. Chodera.
 */

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

}
