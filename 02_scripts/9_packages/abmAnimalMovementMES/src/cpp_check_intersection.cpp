#include <Rcpp.h>
using namespace std;

//' Check if a line segment (step option) crosses a fence
//' @name check_intersection
//' @param origin start location
//' @param target target location
//' @param fence_x1 fence segment x start values
//' @param fence_x2 fence segment x end values
//' @param fence_y1 fence segment y start values
//' @param fence_y2 fence segment y end values
//' @return A Boolean of whether any fence is crossed
//' Useful resource: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/



// CUSTOM STRUCTURES ===========================================================

// declare a Point structure for p1 and q1
struct Point
{
  double x;
  double y;
};

// END STRUCTURES ==============================================================


// HELPER FUNCTIONS ============================================================
bool onSegment(Point p, Point q, Point r) {
  // Given three collinear points p, q, r, the function checks if
  // point q lies on line segment 'pr'
  if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
     q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
    return true;
  return false;
}
int orientation(Point p, Point q, Point r) {
  // To find orientation of ordered triplet (p, q, r).
  // The function returns following values
  // 0 --> p, q and r are collinear
  // 1 --> Clockwise
  // 2 --> Counterclockwise
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  // MAKE SURE THIS IS A DOUBLE AND NOT AN INTEGER!
  double val = ( (q.y - p.y) * (r.x - q.x) ) - ( (q.x - p.x) * (r.y - q.y) );
  if (val == 0) return 0;  // collinear
  return (val > 0)? 1: 2; // clockwise or counterclockwise
}
bool doIntersect(Point p1, Point q1, Point p2, Point q2) {
  // The main function that returns true if line segment 'p1q1'
  // and 'p2q2' intersect.

  // Find the four orientations needed for general and
  // special cases
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);

  // General case
  if (o1 != o2 && o3 != o4) return 1;

  // Special Cases
  // p1, q1 and p2 are collinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1)) return 1;
  // p1, q1 and q2 are collinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1)) return 1;
  // p2, q2 and p1 are collinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2)) return 1;
  // p2, q2 and q1 are collinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2)) return 1;

  return 0; // Doesn't fall in any of the above cases
}

// END HELPER FUNCTIONS ========================================================





// EXPORTED FUNCTION ===========================================================

// [[Rcpp::export]]
double cpp_check_intersection(std::vector<double> origin,
                              std::vector<double> target,
                              std::vector<double> fence_x1,
                              std::vector<double> fence_x2,
                              std::vector<double> fence_y1,
                              std::vector<double> fence_y2) {

  // declare output and within-function data types
  struct Point fence_p, fence_q, p, q;
  bool doesIntersect;
  int L;

  p = {origin[0], origin[1]};
  q = {target[0], target[1]};

  // Loop over fence segments and see if there's a crossing
  L = fence_x1.size();
  for (int i = 0; i < L; i++) {
    fence_p = { fence_x1[i], fence_y1[i] };
    fence_q = { fence_x2[i], fence_y2[i] };
    doesIntersect = doIntersect(p, q, fence_p, fence_q);
    if (doesIntersect) return(true);
  }

  // If there are no crossings
  return(false);
}
// END EXPORTED FUNCTION =======================================================

