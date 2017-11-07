/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.tool;

/**
 *
 * @author sebastien
 */
public class FastTrig {
    
  /** Fast approximation of 1.0 / sqrt(x).
   * See <a href="http://www.beyond3d.com/content/articles/8/">http://www.beyond3d.com/content/articles/8/</a>
   * @param x Positive value to estimate inverse of square root of
   * @return Approximately 1.0 / sqrt(x)
   **/
    public static double invSqrt(double x) {
        final double xhalf = 0.5 * x;
        x = Double.longBitsToDouble(0x5FE6EB50C7B537AAL - (Double.doubleToRawLongBits(x) >> 1));
        return x * (1.5 - xhalf * x * x);
    }

    public static double fast_atan2_unit(double y, double x)  {
        // Normalise such that 0.0 <= y <= x
        final boolean 
                negX = x < 0.0, 
                negY = y < 0.0;
        
        if (negY) {y = -y;}
        if (negX) {x = -x;}
        
        final boolean steep = y > x;
        if (steep) { final double t = x; x = y; y = t; } // Switch

        final double yp = FRAC_BIAS + y;
        final int ind = (int) Double.doubleToRawLongBits(yp);
        
        final double sd = y * COS_TAB[ind] - x * (yp - FRAC_BIAS); // sin(θ-φ) ≡ sinθ cosφ - cosθ sinφ
        
        double θ = ASIN_TAB[ind] + (6.0 + sd * sd) * sd * ONE_SIXTH;
        // Translate back to correct octant
        if (steep) { θ = Math.PI * 0.5 - θ; }
        if (negX)  { θ = Math.PI - θ; }
        if (negY)  { θ = -θ; }

        return θ;
    }

    
  /** Approximation of arctangent.
   *  Slightly faster and substantially less accurate than
   *  {@link Math#atan2(double, double)}.
   **/
  public static double fast_atan2(double y, double x)  {
    final double d2 = x*x + y*y;

    // Bail out if d2 is NaN, zero or subnormal
    if (Double.isNaN(d2) || (Double.doubleToRawLongBits(d2) < 0x10000000000000L)) {
      return Double.NaN;
    }

    // Normalise such that 0.0 <= y <= x
    final boolean negY = y < 0.0;
    if (negY) {y = -y;}
    final boolean negX = x < 0.0;
    if (negX) {x = -x;}
    final boolean steep = y > x;
    if (steep) {
        double t = x;
        x = y;
        y = t;
    }

    // Scale to unit circle (0.0 <= y <= x <= 1.0)
    final double rinv = invSqrt(d2); // rinv ≅ 1.0 / hypot(x, y)
    x *= rinv; // x ≅ cos θ
    y *= rinv; // y ≅ sin θ, hence θ ≅ asin y

    // Hack: we want: ind = floor(y * 256)
    // We deliberately force truncation by adding floating-point numbers whose
    // exponents differ greatly.  The FPU will right-shift y to match exponents,
    // dropping all but the first 9 significant bits, which become the 9 LSBs
    // of the resulting mantissa.
    // Inspired by a similar piece of C code at
    // http://www.shellandslate.com/computermath101.html
    final double yp = FRAC_BIAS + y;
    final int ind = (int) Double.doubleToRawLongBits(yp);

    // Find φ (a first approximation of θ) from the LUT
    final double φ = ASIN_TAB[ind];
    final double cφ = COS_TAB[ind]; // cos(φ)

    // sin(φ) == ind / 256.0
    // Note that sφ is truncated, hence not identical to y.
    final double sφ = yp - FRAC_BIAS;
    final double sd = y * cφ - x * sφ; // sin(θ-φ) ≡ sinθ cosφ - cosθ sinφ

    // asin(sd) ≅ sd + ⅙sd³ (from first 2 terms of Maclaurin series)
    final double d = (6.0 + sd * sd) * sd * ONE_SIXTH;
    
    double θ = φ + d;

    // Translate back to correct octant
    if (steep) { θ = Math.PI * 0.5 - θ; }
    if (negX) { θ = Math.PI - θ; }
    if (negY) { θ = -θ; }

    return θ;
  }

  private static final double ONE_SIXTH = 1.0 / 6.0;
  private static final int FRAC_EXP = 8; // LUT precision == 2 ** -8 == 1/256
  private static final int LUT_SIZE = (1 << FRAC_EXP) + 1;
  private static final double FRAC_BIAS = Double.longBitsToDouble((0x433L - FRAC_EXP) << 52);
  private static final double[] ASIN_TAB = new double[LUT_SIZE];
  private static final double[] COS_TAB = new double[LUT_SIZE];

  static
  {
    /* Populate trig tables */
    for (int ind = 0; ind < LUT_SIZE; ++ ind) {
      double v = ind / (double) (1 << FRAC_EXP);
      double asinv = Math.asin(v);
      COS_TAB[ind] = Math.cos(asinv);
      ASIN_TAB[ind] = asinv;
    }
  }
}
