package com.afieldtrails.common.s2;

import com.google.common.collect.BoundType;
import com.google.common.collect.Range;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2CellUnion;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Utility class encapsulating logic to create {@link Range}s of {@link S2CellId}s.
 */
public class S2RangeUtil {
  private S2RangeUtil() {
  }
  
  /**
   * Return a comparator that can be used for comparing two Range objects of S2CellIds.
   * This method assumes that the ranges being compared will be disjoint and have the same
   * inclusivity on the lower endpoint (closed, specifically).
   */
  public static Comparator<Range<S2CellId>> closedLowerEndpointComparator() {
    return (l, r) -> {
        // Intentional use of == instead of .equals for a faster comparison.
        if (l == r) {
          return 0;
        }
        return l.lowerEndpoint().compareTo(r.lowerEndpoint());
      };
  }

  /**
   * Convert a {@link S2CellUnion} to a list of S2CellId ranges at a maximal level.
   * These ranges are appropriate for id-range based indexing data structures.
   */
  public static void cellIdsToRanges(List<S2CellId> cellIds, int s2Level,
                                     Collection<Range<S2CellId>> outRanges) {
    if (cellIds.size() == 0) {
      return;
    }
    Collections.sort(cellIds);
    S2CellId firstId = cellIds.get(0);
    // If we get ids that are at a more detailed (higher) level than desired, we'll make them
    // more coarse.
    // TODO(lacz): Decide if this is what we actually want. Could be surprising loss of detail.
    if (firstId.level() > s2Level) {
      firstId = firstId.parent(s2Level);
    }
    S2CellId startId = firstId.childBegin(s2Level);
    // endId is after all the children in firstId (it is not contained by firstId)
    S2CellId endId = firstId.childEnd(s2Level);
    for (int i = 1; i < cellIds.size(); ++i) {
      S2CellId nextId = cellIds.get(i);
      // See above about mapping to coarser level.
      if (nextId.level() > s2Level) {
        nextId = nextId.parent(s2Level);
      }
      S2CellId nextIdAtLevel = nextId.childBegin(s2Level);
      if (nextId.id() < endId.id()) {
        // Ignore, we are visiting cells that are at a greater level than the s2 level.
        continue;
      }
      if (endId.equals(nextIdAtLevel)) {
        // Contiguous
        endId = nextId.childEnd(s2Level);
      } else {
        // add range
        outRanges.add(Range.closedOpen(startId, endId));
        startId = nextId.childBegin(s2Level);
        endId = nextId.childEnd(s2Level);
      }
    }
    outRanges.add(Range.closedOpen(startId, endId));
  }

  public static void rangesToCellIds(Iterable<Range<S2CellId>> ranges,
                                     ArrayList<S2CellId> outCellIds) {
    for (Range<S2CellId> range : ranges) {
      rangeToCellIds(range, outCellIds);
    }
  }

  public static void rangeToCellIds(Range<S2CellId> range, ArrayList<S2CellId> outCellIds) {
    S2CellId endId = range.upperEndpoint();
    if (range.upperBoundType().equals(BoundType.OPEN)) {
      endId = endId.prev();
    }
    int endIdLevel = endId.level();
    S2CellId startId = range.lowerEndpoint();
    if (range.lowerBoundType().equals(BoundType.OPEN)) {
      startId = startId.next();
    }
    long remainingLength = endId.id() - startId.id();
    while (remainingLength >= 0) {
      S2CellId cellToAdd = startId;
      // If there's remaining length is long enough, try using a larger cell.
      while (cellToAdd.childPosition(cellToAdd.level()) == 0
          && cellToAdd.parent().childEnd(endIdLevel).lessOrEquals(endId)) {
        cellToAdd = cellToAdd.parent();
      }
      outCellIds.add(cellToAdd);
      startId = cellToAdd.next();
      remainingLength = endId.id() - startId.id();
    }
  }
}
