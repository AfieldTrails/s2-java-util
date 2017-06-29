package com.afieldtrails.common.s2;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S1Interval;
import com.google.common.geometry.S2;
import com.google.common.geometry.S2EdgeUtil;
import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polyline;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;

/**
 * Static utility methods for useful operations on S2Polylines.
 */
public final class S2Polylines {
  private S2Polylines() {
  }

  /**
   * Return whether the points represent a valid S2Polyline.
   */
  public static boolean isValid(@Nonnull List<S2Point> points) {
    // TODO(lacz): do more efficiently.
    return new S2Polyline(new ArrayList<>()).isValid(points);
  }

  /**
   * Construct an {@link S2Polyline} from an ordered collection of
   * {@link S2LatLng}s.
   *
   * @return null if the points do not form a valid polyline.
   */
  @Nullable
  public static S2Polyline fromS2LatLngs(@Nonnull Iterable<S2LatLng> latLngs) {
    ArrayList<S2Point> points = new ArrayList<>();
    for (S2LatLng latLng : latLngs) {
      points.add(latLng.toPoint());
    }
    if (isValid(points)) {
      return new S2Polyline(points);
    } else {
      return null;
    }
  }

  /**
   * Construct a {@link S2Polyline} from an iterable of s2 cell ids.
   */
  public static S2Polyline fromS2CellIds(@Nonnull Iterable<Long> cellIdIterable) {
    ArrayList<S2Point> result = new ArrayList<>();
    for (Long id : cellIdIterable) {
      result.add(S2CellIds.toPoint(id));
    }
    return new S2Polyline(result);
  }

  /**
   * Construct a {@link S2Polyline} from an array of s2 cell ids encoded as longs.
   */
  public static S2Polyline fromS2CellIds(@Nonnull long[] cellIds) {
    ArrayList<S2Point> result = new ArrayList<>();
    for (long id : cellIds) {
      result.add(S2CellIds.toPoint(id));
    }
    return new S2Polyline(result);
  }

  public static S2Polyline subPolyline(S2Polyline polyline, double startAngle, double endAngle) {
    ArrayList<S2Point> points = new ArrayList<>(polyline.numVertices());
    if (startAngle <= 0) {
      points.add(polyline.vertex(0));
    }

    for (int i = 1; i < polyline.numVertices(); ++i) {
      double length = polyline.vertex(i - 1).angle(polyline.vertex(i));
      if (startAngle <= length && startAngle > 0.0) {
        // This code interpolates with respect to arc length rather than
        // straight-line distance, and produces a unit-length result.
        double f = Math.sin(startAngle) / Math.sin(length);
        points.add(S2Point.add(
            S2Point.mul(polyline.vertex(i - 1), (Math.cos(startAngle) - f * Math.cos(length))),
            S2Point.mul(polyline.vertex(i), f)));
      }
      if (!points.isEmpty()) {
        if (endAngle < length) {
          // This code interpolates with respect to arc length rather than
          // straight-line distance, and produces a unit-length result.
          double f = Math.sin(endAngle) / Math.sin(length);
          points.add(S2Point.add(
              S2Point.mul(polyline.vertex(i - 1), (Math.cos(endAngle) - f * Math.cos(length))),
              S2Point.mul(polyline.vertex(i), f)));
          break;
        } else {
          points.add(polyline.vertex(i));
        }
      }
      startAngle -= length;
      endAngle -= length;
    }
    return new S2Polyline(points);
  }

  /**
   * Simplify the input polyline so that the maximum error (l2 distance from a
   * vertex in the input polyline to the resulting polyline) is less than
   * maxError.
   *
   * <p>
   * This algorithm runs in linear (O(nuber of vertices)) time.
   *
   * <p>
   * At each point we compute angles from the source of the current line to
   * projections of the new point at maxError away from the point. A window of
   * angles that are acceptable tolerances for all prior points is kept. If the
   * window doesn't contain the new point, the last seen point is used for the
   * simplified path. Thus, all points are a subset of the points in the input
   * polyline. The first and last points are always retained.
   */
  public static S2Polyline simplify(@Nonnull S2Polyline inputPolyline, @Nonnull S1Angle maxError) {
    // TODO(lacz): Move this to just work on lists of S2Points and make a
    // version for S2Loops.

    int numVerticies = inputPolyline.numVertices();
    ArrayList<S2Point> simplifiedPolyline = new ArrayList<>(numVerticies);
    if (numVerticies <= 2) {
      return inputPolyline;
    }
    S2Point lastEndpoint = inputPolyline.vertex(0);
    simplifiedPolyline.add(lastEndpoint);
    S2Point referencePoint = null;
    S2Point candidateEndpoint = null;
    S1Interval tolerance = S1Interval.full(); // force the first point to be
                                              // added.
    // As we walk away from our last endpoint (A), keep track of the allowable
    // range of angles from
    // that point. We use the first point after the endpoint for a consistent
    // orientation (B).
    // The final point range we much shift orthogonally by maxError to compute
    // the S1Interval for
    // the current point.
    double cosMaxError = Math.cos(maxError.radians());
    double sinMaxError = Math.sin(maxError.radians());
    for (int i = 1; i < numVerticies; ++i) {
      S2Point point = inputPolyline.vertex(i);
      if (point.equals(lastEndpoint)) {
        continue;
      }
      if (referencePoint == null) {
        referencePoint = point;
      }
      if (!tolerance.fastContains(signedAngle(referencePoint, lastEndpoint, point))) {
        simplifiedPolyline.add(candidateEndpoint);
        lastEndpoint = candidateEndpoint;
        tolerance = S1Interval.full();
        referencePoint = point;
        candidateEndpoint = null;
      }
      S2Point directVector = S2Point.minus(point, lastEndpoint);
      S2Point scaledOrtho = S2Point.mul(S2Point.normalize(S2Point.crossProd(point, directVector)),
          sinMaxError);
      S2Point scaledPoint = S2Point.mul(point, cosMaxError);
      S2Point offset0 = S2Point.add(scaledPoint, scaledOrtho);
      S2Point offset1 = S2Point.add(scaledPoint, S2Point.neg(scaledOrtho));
      S1Interval pointTolerance = S1Interval.fromPointPair(
          signedAngle(referencePoint, lastEndpoint, offset0),
          signedAngle(referencePoint, lastEndpoint, offset1));
      tolerance = tolerance.intersection(pointTolerance);
      candidateEndpoint = point;
    }
    simplifiedPolyline.add(inputPolyline.vertex(numVerticies - 1));
    return new S2Polyline(simplifiedPolyline);
  }

  /**
   * Returns the 'signed' interior angle of B for the angle ABC.
   * 
   * The result is positive if the points are counter-clockwise and negative if
   * the points are clockwise.
   */
  private static double signedAngle(@Nonnull S2Point a, @Nonnull S2Point b, @Nonnull S2Point c) {
    double angle = S2.angle(a, b, c);
    return (S2.robustCCW(a, b, c) > 0) ? angle : -angle;
  }

  /**
   * Builds and initializes an {@link S2EdgeUtil.XYZPruner} for quickly testing
   * possible intersections with the input polyline.
   */
  public static S2EdgeUtil.XYZPruner buildXYZPruner(@Nonnull S2Polyline polyline) {
    S2EdgeUtil.XYZPruner pruner = new S2EdgeUtil.XYZPruner();
    for (int i = 0; i < polyline.numVertices() - 1; ++i) {
      pruner.addEdgeToBounds(polyline.vertex(i), polyline.vertex(i + 1));
    }
    return pruner;
  }

  /**
   * Gives the S2Points for all intersections between polylineA and polylineB.
   * Returns
   */
  public static HashSet<S2Point> findIntersections(@Nonnull S2Polyline polylineA,
                                                   @Nonnull S2Polyline polylineB) {
    // TODO(lacz): This should work, but an edge index might be more efficient.
    HashSet<S2Point> outIntersections = new HashSet<>();
    S2Polyline shorterPolyline = (polylineA.numVertices() < polylineB.numVertices())
        ? polylineA : polylineB;
    S2Polyline longerPolyline = (polylineA.numVertices() < polylineB.numVertices())
        ? polylineB : polylineA;

    S2EdgeUtil.XYZPruner pruner = buildXYZPruner(shorterPolyline);
    ArrayList<Integer> candidates = new ArrayList<>();
    pruner.setFirstIntersectPoint(longerPolyline.vertex(0));
    for (int i = 1; i < longerPolyline.numVertices(); ++i) {
      if (pruner.intersects(longerPolyline.vertex(i))) {
        candidates.add(i - 1);
      }
    }
    for (Integer edgeIndex : candidates) {
      Preconditions.checkState(longerPolyline.numVertices() > edgeIndex + 1);
      S2Point a0 = longerPolyline.vertex(edgeIndex);
      S2Point a1 = longerPolyline.vertex(edgeIndex + 1);
      S2Point b0 = shorterPolyline.vertex(0);
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a0, a1,
          shorterPolyline.vertex(0));
      for (int j = 1; j < shorterPolyline.numVertices(); ++j) {
        S2Point b1 = shorterPolyline.vertex(j);
        switch (crosser.robustCrossing(b1)) {
          case -1:
            // No intersection.
            break;
          case 0:
            // Intersection at an endpoint.
            // Only use one of the endpoints (a0 in this case) to avoid double additions.
            if (a0.equals(b0) || a0.equals(b1)) {
              outIntersections.add(a0);
            }
            if (a1.equals(b0) || a1.equals(b1)) {
              outIntersections.add(a1);
            }
            break;
          case 1:
            // A robust intersection.
            // There is an intersection between a and b.
            outIntersections.add(S2EdgeUtil.getIntersection(a0, a1, b0, b1));
            break;
        }
        b0 = b1;
      }
    }
    return outIntersections;
  }

  /**
   * Return the arc length (in radians) of the polyline from the startIndex to the endIndex
   * (inclusive).
   */
  public static double getArcLengthOfSubpolyline(@Nonnull S2Polyline polyline,
                                                 int startIndex, int endIndex) {
    Preconditions.checkArgument(startIndex <= endIndex);
    Preconditions.checkArgument(startIndex >= 0);
    Preconditions.checkArgument(endIndex < polyline.numVertices());
    double radians = 0.0;
    for (int i = startIndex; i < endIndex; ++i) {
      radians += polyline.vertex(i).angle(polyline.vertex(i + 1));
    }
    return radians;
  }

  /**
   * Return the arc length (in radians) from the start of a polyline to the given point.  The point
   * must be on the polyline or this value will be meaningless.
   */
  public static double getArcLengthToPoint(@Nonnull S2Polyline polyline, @Nonnull S2Point point) {
    int edgeIndex = polyline.getNearestEdgeIndex(point);
    double radians = getArcLengthOfSubpolyline(polyline, 0, edgeIndex);
    S2Point onPolyline = S2EdgeUtil.getClosestPoint(point, polyline.vertex(edgeIndex),
        polyline.vertex(edgeIndex + 1));
    radians += polyline.vertex(edgeIndex).angle(onPolyline);
    return radians;
  }

  public static S2Point getClosestPoint(@Nonnull S2Polyline polyline, @Nonnull S2Point point) {
    // TODO(lacz): Use S2Polylines distance to point (slightly more efficient)
    int nearestEdgeIndex = polyline.getNearestEdgeIndex(point);
    return S2EdgeUtil.getClosestPoint(point, polyline.vertex(nearestEdgeIndex),
        polyline.vertex(nearestEdgeIndex + 1));
  }
  
  /**
   * Return the distance (in radians) between the given polyline and the point.
   */
  public static double getDistanceToPoint(@Nonnull S2Polyline polyline, @Nonnull S2Point point) {
    Preconditions.checkState(polyline.numVertices() > 0, "Empty polyline");

    if (polyline.numVertices() == 1) {
      return polyline.vertex(0).angle(point);
    }

    // Initial value larger than any possible distance on the unit sphere.
    S1Angle minDistance = S1Angle.radians(10);

    // Find the line segment in the polyline that is closest to the point given.
    for (int i = 0; i < polyline.numVertices() - 1; ++i) {
      S1Angle distanceToSegment = S2EdgeUtil.getDistance(point, polyline.vertex(i),
          polyline.vertex(i + 1));
      if (distanceToSegment.lessThan(minDistance)) {
        minDistance = distanceToSegment;
      }
    }
    return minDistance.radians();
  }
}
