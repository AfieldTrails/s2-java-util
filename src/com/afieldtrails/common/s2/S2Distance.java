package com.afieldtrails.common.s2;

import com.google.common.geometry.S2LatLng;

/**
 * Utility class for working with distances and S2 geometry.
 */
public class S2Distance {
  private S2Distance() {}

  private static final double WGS84_MAJOR_AXIS = 6378137.0;
  private static final double WGS84_SEMI_MAJOR_AXIS = 6356752.3142;

  // The radius assumed by s2 geometry.
  // Should match S2LatLng.EARTH_RADIUS_METERS
  private static final double S2_EARTH_RADIUS_METERS = 6367000.0;


  public static double metersToRadians(double meters) {
    return meters / S2_EARTH_RADIUS_METERS;
  }

  /**
   * Compute the approximate distance between two points in the WGS84 ellipsoid.
   */
  public static float computeWGS84DistanceMeters(S2LatLng latLng1, S2LatLng latLng2) {
    return computeWGS84DistanceMeters(latLng1.latRadians(),
        latLng1.lngRadians(),
        latLng2.latRadians(), latLng2.lngRadians());
  }

  /**
   * Compute the approximate distance between two points in the WGS84 ellipsoid.
   * Latitude and Longitude values are in radians.
   */
  private static float computeWGS84DistanceMeters(double lat1,
        double lon1, double lat2, double lon2) {
    // Copied from Location.getDistance so that it can be executed in unit tests (and avoid a
    // couple superfluous conversions)
    // Based on http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
    // using the "Inverse Formula" (section 4)

    int MAXITERS = 20;

    double a = WGS84_MAJOR_AXIS;
    double b = WGS84_SEMI_MAJOR_AXIS;
    double f = (a - b) / a;
    double aSqMinusBSqOverBSq = (a * a - b * b) / (b * b);

    double L = lon2 - lon1;
    double A = 0.0;
    double U1 = Math.atan((1.0 - f) * Math.tan(lat1));
    double U2 = Math.atan((1.0 - f) * Math.tan(lat2));

    double cosU1 = Math.cos(U1);
    double cosU2 = Math.cos(U2);
    double sinU1 = Math.sin(U1);
    double sinU2 = Math.sin(U2);
    double cosU1cosU2 = cosU1 * cosU2;
    double sinU1sinU2 = sinU1 * sinU2;

    double sigma = 0.0;
    double deltaSigma = 0.0;
    double cosSqAlpha = 0.0;
    double cos2SM = 0.0;
    double cosSigma = 0.0;
    double sinSigma = 0.0;
    double cosLambda = 0.0;
    double sinLambda = 0.0;

    double lambda = L; // initial guess
    for (int iter = 0; iter < MAXITERS; iter++) {
      double lambdaOrig = lambda;
      cosLambda = Math.cos(lambda);
      sinLambda = Math.sin(lambda);
      double t1 = cosU2 * sinLambda;
      double t2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
      double sinSqSigma = t1 * t1 + t2 * t2; // (14)
      sinSigma = Math.sqrt(sinSqSigma);
      cosSigma = sinU1sinU2 + cosU1cosU2 * cosLambda; // (15)
      sigma = Math.atan2(sinSigma, cosSigma); // (16)
      double sinAlpha = (sinSigma == 0) ? 0.0 :
                        cosU1cosU2 * sinLambda / sinSigma; // (17)
      cosSqAlpha = 1.0 - sinAlpha * sinAlpha;
      cos2SM = (cosSqAlpha == 0) ? 0.0 :
               cosSigma - 2.0 * sinU1sinU2 / cosSqAlpha; // (18)

      double uSquared = cosSqAlpha * aSqMinusBSqOverBSq; // defn
      A = 1 + (uSquared / 16384.0) * // (3)
          (4096.0 + uSquared *
              (-768 + uSquared * (320.0 - 175.0 * uSquared)));
      double B = (uSquared / 1024.0) * // (4)
          (256.0 + uSquared *
              (-128.0 + uSquared * (74.0 - 47.0 * uSquared)));
      double C = (f / 16.0) *
          cosSqAlpha *
          (4.0 + f * (4.0 - 3.0 * cosSqAlpha)); // (10)
      double cos2SMSq = cos2SM * cos2SM;
      deltaSigma = B * sinSigma * // (6)
          (cos2SM + (B / 4.0) *
              (cosSigma * (-1.0 + 2.0 * cos2SMSq) -
                  (B / 6.0) * cos2SM *
                      (-3.0 + 4.0 * sinSigma * sinSigma) *
                      (-3.0 + 4.0 * cos2SMSq)));

      lambda = L +
          (1.0 - C) * f * sinAlpha *
              (sigma + C * sinSigma *
                  (cos2SM + C * cosSigma *
                      (-1.0 + 2.0 * cos2SM * cos2SM))); // (11)

      double delta = (lambda - lambdaOrig) / lambda;
      if (Math.abs(delta) < 1.0e-12) {
        break;
      }
    }

    float distance = (float) (b * A * (sigma - deltaSigma));
    return distance;
    /*
    // TODO(lacz): Include a version to return initial and final bearings.
    results.mDistance = distance;
    float initialBearing = (float) Math.atan2(cosU2 * sinLambda,
        cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);
    initialBearing *= 180.0 / Math.PI;
    results.mInitialBearing = initialBearing;
    float finalBearing = (float) Math.atan2(cosU1 * sinLambda,
        -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);
    finalBearing *= 180.0 / Math.PI;
    results.mFinalBearing = finalBearing;
    results.mLat1 = lat1;
    results.mLat2 = lat2;
    results.mLon1 = lon1;
    results.mLon2 = lon2;
    */
  }
}
