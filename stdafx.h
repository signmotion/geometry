#pragma once

// @error Если при компиляции возникает ошибка >
//        http://stackoverflow.com/questions/3376224/ms-vc-iostream-compile-error

#define _USE_MATH_DEFINES

#include <iostream>
//#include <deque>
#include <vector>
//#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp> 
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/concepts/check.hpp>
#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp>
//#include <boost/tuple/tuple_io.hpp>


using namespace std;

namespace bg = boost::geometry;
namespace bgm = bg::model;
