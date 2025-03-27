#pragma once

#include "base/Global.h"
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

    // Builds the RTree
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

class SpatialDataReader 
{
  public:
    /**     **/
    struct Pt 
    {
      double x, y, z, val;
      Pt(double x_, double y_, double z_, double val_=0);
      double dist(const Pt& m);
      void print() const;
    };
    /**     **/

    SpatialDataReader( string fn_ );
    const vector<Pt>& get_points() const { return points; };
    void printPoints() const;
    double value_at( double x, double y, double z, uint n_nb = 2 );

  private:
    void readFromFile();
    
    vector<Pt> points;
    string fn;

    using RTreePoint = bg::model::point<double, 3, bg::cs::cartesian>;
    bgi::rtree<pair<RTreePoint, Pt>, bgi::quadratic<16>> rtree;

    // Burocratic stuff
    friend ostream& operator<<(ostream& os, const SpatialDataReader & m);
    friend ostream& operator<<(ostream& os, const SpatialDataReader::Pt & m);
};

ostream& operator<<(ostream& os, const SpatialDataReader & m);
ostream& operator<<(ostream& os, const SpatialDataReader::Pt & m);
