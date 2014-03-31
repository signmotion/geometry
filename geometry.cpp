#include "stdafx.h"
#include "external/triangle/triangle.h"


/**
* ������ �������, �.
*/
const double WORLD_PLANET_RADIUS = 5000.0 * 1000.0;

/**
* ������ - ���������� �������, � ������� �������� �������������.
* ������� � ������.
*/
const double SECTOR_SIZE = 1000.0;
const double SECTOR_HALF_SIZE = SECTOR_SIZE / 2.0;


//typedef std::pair< double, double>  coord_xy_t;
typedef bgm::d2::point_xy< double >  p_t;
typedef bgm::polygon< p_t >          polygon_t;
typedef boost::tuple< p_t, p_t, p_t >  triangle_t;



/**
* @return ������� �����. ������������, �� �����������, ������ �������.
*/
polygon_t drawCircle( double radius, size_t vertex ) {
	
	assert( (radius >= 100.00) && "����� ��������� �� ����, �������� >= 100.00" );

    polygon_t circle;
    const auto a = M_PI / vertex * 2;
    for (double theta = 0; theta <= (2.0 * M_PI); theta += a) {
		const auto x = radius * cos( theta );
		const auto y = radius * sin( theta );
        const auto c = bg::make< p_t >( floor( x ), floor( y ) );
        bg::append( circle, c );
    }
	bg::correct( circle );

    return circle;
}




/**
* ����������� ����� � ���������������. ���� ������� ������ ��������������.
*/
void testIntersectionCircleBox() {
    
    // �������
	// (!) ������� ����� �������� ��� ����� ����� �������� � ����, ���
	//     ������ (��. ����) � ���� �� �����������.
    const bgm::polygon< p_t >  circle = drawCircle( 5e6, 20 );
    cout << "������� ������� " << bg::dsv( circle ) << endl;
    cout << "������� " << bg::area( circle ) << endl << endl;


    // ������
    const auto aimCoord = p_t( 0.0, WORLD_PLANET_RADIUS );
    bgm::box< p_t >  sector(
        p_t( aimCoord.x() - SECTOR_HALF_SIZE,  aimCoord.y() - SECTOR_HALF_SIZE ),
        p_t( aimCoord.x() + SECTOR_HALF_SIZE,  aimCoord.y() + SECTOR_HALF_SIZE )
    );
    bg::correct( sector );
    cout << "������� ������� " << bg::dsv( sector ) << endl;
    cout << "������� " << bg::area( sector ) << endl << endl;


    // �����������
    std::vector< polygon_t > pv;
    bg::intersection( circle, sector, pv );
    cout << "���������� ��������� ����� ����������� " << pv.size() << endl << endl;
	assert( (pv.size() == 1) && "�������� 1 �������." );

    // ����������� �������
    for (auto itr = pv.cbegin(); itr != pv.cend(); ++itr) {
        const polygon_t poly = *itr;
        // ������ ������� - ��� ��������� �������� ��������
        cout << "������� ������� � ������� " << bg::dsv( poly ) << endl;
        cout << "������� " << bg::area( poly ) << endl << endl;
        // ������� ��������
        const auto vertex = poly.outer();
        for (auto itrVertex = vertex.cbegin(); itrVertex != vertex.cend(); ++itrVertex) {
            const auto t = *itrVertex;
            cout << t.x() << " " << t.y() << endl;
        }
        cout << endl;
    }

}







/**
* @helper testBurnTrianglePolygonAndTriangulateResult()
*/
void printReportTriangulate(
    struct triangulateio* io,
    int markers,
    int reporttriangles,
    int reportneighbors,
    int reportsegments,
    int reportedges,
    int reportnorms
) {
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist ? io->pointlist[i * 2 + j] : 0);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g", io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist ? io->pointmarkerlist[i] : 0);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist ? io->trianglelist[i * io->numberofcorners + j] : 0);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist ? io->triangleattributelist[i * io->numberoftriangleattributes + j] : 0);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist ? io->neighborlist[i * 3 + j] : 0);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist ? io->segmentlist[i * 2 + j] : 0);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist ? io->segmentmarkerlist[i] : 0);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist ? io->edgelist[i * 2 + j] : 0);
      }
      if (reportnorms && io->edgelist && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist ? io->normlist[i * 2 + j] : 0);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist ? io->edgemarkerlist[i] : 0);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}






/**
* ��������� ������� �� ������������. ������� ����� ��������� ����.
*
* @source http://www.cs.cmu.edu/~quake/triangle.html 
* @source external/triangle/tricall.c
* @source http://people.sc.fsu.edu/~jburkardt/c_src/triangle/triangle.html
* @see �������� ���������� ����, ����� � triangulate()
*
* ������ ����������� ��� ����� ���� �����. ������ - ��� ����:
* 1. ������ ������������ ����� ������� ������ (��� ���).
*    � ���� ������, ��������� ������ ������ ������ 'pointlist'.
*
* 2. ������ - ������ � ������.
*    ������ 'pointlist' �����������:
*      �) ������ ��������� 'segmentlist'; ������������ ����� ������ ������,
*         ����������� ��� ������ (������ ���� ���������� �������)
*         @example http://people.sc.fsu.edu/~jburkardt/data/triangle_files/double_hex2.poly
*                  (� ������� ������� ������ ������ 8-��� ������ ��� ������ �������)
*      �) ����; ������������� ��� ���������� ������ �������� � 'segmentlist' �����
*
* @helper testBurnTrianglePolygonAndTriangulateResult()
*/
std::vector< triangle_t > triangulatePolygon( polygon_t& poly ) {

    // ������� ������ ��������� ������������ ������ �������
    bg::correct( poly );

    // ���������� ���������, ����������� ��� ������ ������ ������������
    struct triangulateio in;
    struct triangulateio out;

    cout << "������� ��� ������������, ������� ������: " <<
        poly.outer().size() << " �������" << endl << endl;
    cout << "������� ��� ������������, ���������� ������: " <<
        poly.inners().size() << endl << endl;

	// ������� ������, ������ ����
    const auto outerVertex = poly.outer();
	// ��������� ������� ��������� ������, -1
	const auto countOuterVertex = outerVertex.size() - 1;

	// ���������� ������ (����), ����� ������������� ��� ���� ���������
    const auto innerPoly = poly.inners();
	int countInnerVertex = 0;
	for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
        const auto hole = *itrPoly;
		// ��������� ������� ��������� ������, -1
		countInnerVertex += (int)hole.size() - 1;
	}


	const auto countVertex = countOuterVertex + countInnerVertex;

	cout << endl << "��������� 'pointlist'" << endl;

    in.numberofpoints = (int)countVertex;
    // ������� ������� ����� ���������
    in.pointlist = new real_t[ in.numberofpoints * 2 ];
	in.numberofpointattributes = 0;
    in.pointattributelist = nullptr;
    in.pointmarkerlist = nullptr;

    int pointerPoint2 = 0;

	// ������� ������ (��������� ������� ������ ��������� ������, -1) ...
    for (auto itrVertex = outerVertex.cbegin(); itrVertex != outerVertex.cend() - 1; ++itrVertex) {
        const auto vertex = *itrVertex;
        const auto x = (real_t)vertex.x();
        const auto z = (real_t)vertex.y();
        cout << "������� ������, ������� " << x << ", " << z << endl;
        in.pointlist[pointerPoint2 + 0] = x;
        in.pointlist[pointerPoint2 + 1] = z;
        pointerPoint2 += 2;
    }

	// ... ���� ���������� ������
	for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
        const auto hole = *itrPoly;
		// (��������� ������� ���� ��������� ������, -1)
		for (auto itrVertex = hole.cbegin(); itrVertex != hole.cend() - 1; ++itrVertex) {
			const auto vertex = *itrVertex;
			const auto x = (real_t)vertex.x();
			const auto z = (real_t)vertex.y();
			cout << "���������� ������, ������� " << x << ", " << z << endl;
			in.pointlist[pointerPoint2 + 0] = x;
			in.pointlist[pointerPoint2 + 1] = z;
			pointerPoint2 += 2;
		}

	} // for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly)


	// �������� � ���� ��������� ������ ��� ������� ���
	if (countInnerVertex == 0) {

		in.numberofsegments = in.numberofholes = 0;
		in.segmentlist = in.segmentmarkerlist = nullptr;
		in.holelist = nullptr;

	} else {

		cout << endl << "��������� 'segmentlist'" << endl;

		// ������ ������ ����������
		in.numberofsegments = (int)countVertex;
		in.segmentlist = new int[ in.numberofsegments * 2 ];
		in.segmentmarkerlist = nullptr;

		int pointerSegment2 = 0;
		int numberPoint = 0;
		int firstNumberPoint = numberPoint;

		// ������� ������ � ���������� ������� ...
		cout << endl << "������� ������" << endl;
		for (auto itrVertex = outerVertex.cbegin(); itrVertex != outerVertex.cend() - 2; ++itrVertex) {
			in.segmentlist[pointerSegment2 + 0] = numberPoint;
			in.segmentlist[pointerSegment2 + 1] = numberPoint + 1;
			cout << "������� ������, ������� " <<
				in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
			++numberPoint;
			pointerSegment2 += 2;
		}
		in.segmentlist[pointerSegment2 + 0] = numberPoint;
		in.segmentlist[pointerSegment2 + 1] = firstNumberPoint;
		cout << "������� ������, ���������� ������� " <<
			in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
		++numberPoint;
		pointerSegment2 += 2;

		// ... ���� ���������� ������ � ����������� ���������
		for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
			const auto hole = *itrPoly;
			cout << endl << "���������� ������" << endl;
			firstNumberPoint = numberPoint;
			for (auto itrVertex = hole.cbegin(); itrVertex != hole.cend() - 2; ++itrVertex) {
				in.segmentlist[pointerSegment2 + 0] = numberPoint;
				in.segmentlist[pointerSegment2 + 1] = numberPoint + 1;
				cout << "���������� ������, ������� " <<
					in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
				++numberPoint;
				pointerSegment2 += 2;
			}
			in.segmentlist[pointerSegment2 + 0] = numberPoint;
			in.segmentlist[pointerSegment2 + 1] = firstNumberPoint;
			cout << "���������� ������, ���������� ������� " <<
				in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
			++numberPoint;
			pointerSegment2 += 2;

		} // for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly)


		// ��������� ������������, ����� �������� �������� ������ (holes)
		in.numberofholes = (int)innerPoly.size();
		in.holelist = new real_t[ in.numberofholes * 2 ];

		int pointerHole2 = 0;

		for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
			const auto hole = *itrPoly;
			p_t center;
			bg::centroid( hole, center );
			cout << endl << "����, ����� " << center.x() << ", " << center.y() << endl << endl;
			in.holelist[pointerHole2 + 0] = center.x();
			in.holelist[pointerHole2 + 1] = center.y();
			pointerHole2 += 2;
		}

	} // else if (countInnerVertex == 0)


	// ������
	in.numberofregions = 0;
	in.regionlist = nullptr;


    printf( "Input point set:\n\n" );
    printReportTriangulate( &in, 0, 0, 0, 1, 0, 0 );


    // �������������� ���������, � ������� ����� ����� ���������
    out.pointlist = 
        out.pointattributelist = 
		out.triangleattributelist = 
		out.trianglearealist = 
	nullptr;
	out.pointmarkerlist = 
		out.trianglelist = 
		out.neighborlist = 
		out.segmentlist = 
		out.segmentmarkerlist = 
		out.edgelist = 
		out.edgemarkerlist =
	nullptr;

	/**
	* @see params > http://www.cs.cmu.edu/~quake/triangle.switch.html
		-p* Triangulates a Planar Straight Line Graph (.poly file).
		-r* Refines a previously generated mesh.
		-q* Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.
		-a Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file.
		-u Imposes a user-defined constraint on triangle size.
		-A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
		-c* Encloses the convex hull with segments.
		-D* Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.
		-j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
		-e Outputs (to an .edge file) a list of edges of the triangulation.
		-v Outputs the Voronoi diagram associated with the triangulation. Does not attempt to detect degeneracies, so some Voronoi vertices may be duplicated.
		-n Outputs (to a .neigh file) a list of triangles neighboring each triangle.
		-g Outputs the mesh to an Object File Format (.off) file, suitable for viewing with the Geometry Center's Geomview package.
		-B Suppresses boundary markers in the output .node, .poly, and .edge output files.
		-P Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain constraining segments on later refinements of the mesh.
		-N Suppresses the output .node file.
		-E Suppresses the output .ele file.
		-I Suppresses mesh iteration numbers.
		-O Suppresses holes: ignores the holes in the .poly file.
		-X* Suppresses exact arithmetic.
		-z* Numbers all items starting from zero (rather than one). Note that this switch is normally overrided by the value used to number the first vertex of the input .node or .poly file. However, this switch is useful when calling Triangle from another program.
		-o2 Generates second-order subparametric elements with six nodes each.
		-Y Prohibits the insertion of Steiner points on the mesh boundary. If specified twice (-YY), it prohibits the insertion of Steiner points on any segment, including internal segments.
		-S Specifies the maximum number of added Steiner points.
		-i* Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
		-F* Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
		-l Uses only vertical cuts in the divide-and-conquer algorithm. By default, Triangle uses alternating vertical and horizontal cuts, which usually improve the speed except with vertex sets that are small or short and wide. This switch is primarily of theoretical interest.
		-s Specifies that segments should be forced into the triangulation by recursively splitting them at their midpoints, rather than by generating a constrained Delaunay triangulation. Segment splitting is true to Ruppert's original algorithm, but can create needlessly small triangles. This switch is primarily of theoretical interest.
		-C* Check the consistency of the final mesh. Uses exact arithmetic for checking, even if the -X switch is used. Useful if you suspect Triangle is buggy.
		-Q* Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.
		-V* Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail. `-V' gives information on algorithmic progress and detailed statistics.
		-h Help: Displays complete instructions.
	  * ������� �������� ��������/���������� ���������.
	*/
	// @optimize '-X'
	// @optimize '-Q'
	// @optimize ��� �������: '-i' ��� '-F'?
	// @interesting '-V'
    triangulate( "Xiz", &in, &out, nullptr );

    printf( "Initial triangulation:\n\n" );
    printReportTriangulate( &out, 0, 1, 0, 1, 0, 0 );

    // �������� ������������ ������������
	// (!) triangle() �� ��������� �� �������� '-c' � ������ ������
	// �� ������ �������� �������������. �.�. ����� ��� ����� �����
	// �������������, ������� ����������� �������� ������, ��������
	// � ����� �� �������������� ��������.
    assert( (out.trianglelist && (out.numberoftriangles > 0) ) && "������ �� ������� �� ������������." );
    assert( (out.numberofcorners == 3) && "������� ������ ������������." );
    std::vector< triangle_t >  vt;
    for (int i = 0; i < out.numberoftriangles; ++i) {
        const auto q = i * out.numberofcorners;

        // �������� ����� �����
        auto pn = out.trianglelist[ q + 0 ];
        // �������� ���������� �����
        auto x = out.pointlist[ pn * 2 + 0 ];
        auto z = out.pointlist[ pn * 2 + 1 ];
		const p_t coord1 = bg::make< p_t >( x, z );

        pn = out.trianglelist[ q + 1 ];
        x = out.pointlist[ pn * 2 + 0 ];
        z = out.pointlist[ pn * 2 + 1 ];
		const p_t coord2 = bg::make< p_t >( x, z );

        pn = out.trianglelist[ q + 2 ];
        x = out.pointlist[ pn * 2 + 0 ];
        z = out.pointlist[ pn * 2 + 1 ];
		const p_t coord3 = bg::make< p_t >( x, z );

		// ���������, ��� ����������� ����������� �������� ������
		// �������� ������ ����� ����� ���� ������������
		const auto c = bg::make< p_t >(
			(coord1.x() + coord2.x() + coord3.x()) / 3.0,
		    (coord1.y() + coord2.y() + coord3.y()) / 3.0
	    );
		if ( bg::within( c, poly ) ) {
			const triangle_t tri( coord1, coord2, coord3 );
			vt.push_back( tri );
		} else {
			cout << "����������� " <<
				bg::dsv( coord1 ) << " " <<
				bg::dsv( coord2 ) << " " <<
				bg::dsv( coord3 ) << " " << " ��������" <<
			endl << endl;
		}

    } // for (int i = 0; i < out.numberoftriangles; ++i)

    cout << "������� ������������� " << vt.size() << endl;


    // ������� �� �����
    delete[] in.pointlist;
    delete[] in.pointattributelist;
    delete[] in.pointmarkerlist;
    delete[] in.regionlist;
    delete[] out.pointlist;
    delete[] out.pointattributelist;
    delete[] out.pointmarkerlist;
    delete[] out.trianglelist;
    delete[] out.triangleattributelist;
    delete[] out.trianglearealist;
    delete[] out.neighborlist;
    delete[] out.segmentlist;
    delete[] out.segmentmarkerlist;
    delete[] out.edgelist;
    delete[] out.edgemarkerlist;

    return vt;
}








/**
* ����������� ������������ � ������ (������� ��� �����). ������ �����
* ������ ������������ � ����� ������ ������������, �������� � �� ����.
* ����� ���������� ������ ����������� ������ � ������������.
*/
void testBurnTrianglePolygonAndTriangulateResult1() {

	// �����������
	/*
	triangle_t tri(
		bg::make< p_t >( 0.0,     0.0 ),
		bg::make< p_t >( 0.0,   100.0 ),
		bg::make< p_t >( 100.0,   0.0 )
	);
	*/
	polygon_t tri;
	bg::append( tri, bg::make< p_t >(   0.0,   0.0 ) );
	bg::append( tri, bg::make< p_t >(   0.0, 100.0 ) );
	bg::append( tri, bg::make< p_t >( 100.0,   0.0 ) );
	bg::correct( tri );
	cout << "�����������: " << bg::dsv( tri ) << endl;

	// ������ ������ ������������
	polygon_t ring;
	bg::append( ring, bg::make< p_t >( 10.0, 10.0 ) );
	bg::append( ring, bg::make< p_t >( 10.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 10.0 ) );
	bg::correct( ring );
	cout << "������ ������ ������������ (����): " << bg::dsv( ring ) << endl;


	cout << endl << "�������� � ������������ ���� �� ����� ������." << endl << endl;
	// @see http://www.boost.org/doc/libs/1_47_0/libs/geometry/doc/html/geometry/reference/algorithms/difference.html
    std::vector< polygon_t > figure;
    bg::difference( tri, ring, figure );
	
	cout << endl << "��������� ����� ��������� ����: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "������� " << bg::dsv( *itr ) << endl;
    }

	assert( (figure.size() == 1) && "������ ���������� 1 �������." );

	// ����� ���������� ������ �� ������������
	polygon_t firstPoly = *figure.cbegin();
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "������ ���������� ����� �������������." );

}







/**
* @see testBurnTrianglePolygonAndTriangulateResult1()
*
* � ������������ �������� ��� ����������������� ����.
*/
void testBurnTrianglePolygonAndTriangulateResult2() {

	// �����������
	polygon_t tri;
	bg::append( tri, bg::make< p_t >(   0.0,   0.0 ) );
	bg::append( tri, bg::make< p_t >(   0.0, 100.0 ) );
	bg::append( tri, bg::make< p_t >( 100.0,   0.0 ) );
	bg::correct( tri );
	cout << "�����������: " << bg::dsv( tri ) << endl;

	// ������ 1 ������ ������������
	polygon_t ring1;
	bg::append( ring1, bg::make< p_t >( 10.0, 10.0 ) );
	bg::append( ring1, bg::make< p_t >( 10.0, 50.0 ) );
	bg::append( ring1, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring1, bg::make< p_t >( 20.0, 10.0 ) );
	bg::correct( ring1 );
	cout << "������ 1 ������ ������������ (����): " << bg::dsv( ring1 ) << endl;

	// ������ 2 ������ ������������, �� ����������� ������ 1
	polygon_t ring2;
	bg::append( ring2, bg::make< p_t >( 1.0, 1.0 ) );
	bg::append( ring2, bg::make< p_t >( 1.0, 5.0 ) );
	bg::append( ring2, bg::make< p_t >( 2.0, 5.0 ) );
	bg::append( ring2, bg::make< p_t >( 2.0, 2.0 ) );
	bg::correct( ring2 );
	cout << "������ 2 ������ ������������ (����): " << bg::dsv( ring2 ) << endl;
	assert( !bg::intersects( ring1, ring2) && "������ �� ������ ������������." );

	cout << endl << "�������� � ������������ ���� �� ����� ������ 1." << endl;
    std::vector< polygon_t > figure;
    bg::difference( tri, ring1, figure );
	assert( (figure.size() == 1) && "������ ���������� 1 �������." );

	cout << endl << "�������� � ������������ ���� �� ����� ������ 2." << endl;
	polygon_t firstPoly = *figure.cbegin();
	figure.clear();
	bg::difference( firstPoly, ring2, figure );
	assert( (figure.size() == 1) && "������ ���������� 1 �������." );
	
	cout << endl << "��������� ����� ��������� 2 ���: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "������� " << bg::dsv( *itr ) << endl;
    }

	// ����� ���������� ������ �� ������������
	firstPoly = *figure.cbegin();
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "������ ���������� ����� �������������." );

}







/**
* @see testBurnTrianglePolygonAndTriangulateResult1()
*
* � ������������ �������� ����, �������� �� ���� �������� �������.
*/
void testBurnTrianglePolygonAndTriangulateResult3() {

	// �����������
	polygon_t tri;
	bg::append( tri, bg::make< p_t >(   0.0,   0.0 ) );
	bg::append( tri, bg::make< p_t >(   0.0, 100.0 ) );
	bg::append( tri, bg::make< p_t >( 100.0,   0.0 ) );
	bg::correct( tri );
	cout << "�����������: " << bg::dsv( tri ) << endl;

	// ������, ����� - ������ ������������, ����� - �������
	polygon_t ring;
	bg::append( ring, bg::make< p_t >( -5.0, -5.0 ) );
	bg::append( ring, bg::make< p_t >( -5.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, -5.0 ) );
	bg::correct( ring );
	cout << "������ �������� ����� ������������: " << bg::dsv( ring ) << endl;
	//assert( bg::overlaps( tri, ring ) && "����������� � ������ ������ �������� ����������� ���� �����." );
	//assert( bg::intersects( tri, ring ) && !bg::within( ring, tri ) && "����������� � ������ ������ �������� ����������� ���� �����." );


	cout << endl << "�������� � ������������ ����� �� ����� ������." << endl << endl;
    std::vector< polygon_t > figure;
    bg::difference( tri, ring, figure );
	
	cout << endl << "��������� ����� ���������: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "������� " << bg::dsv( *itr ) << endl;
    }

	assert( (figure.size() == 1) && "������ ���������� 1 �������." );

	// ����� �������� ������� �� ������������
	polygon_t firstPoly = *figure.cbegin();
	assert( (firstPoly.inners().size() == 0) && "������ ������ ���� � ��������, �� ��� ����." );
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "������ ���������� ����� �������������." );

}







int main() {

	setlocale( LC_ALL, "Russian" );
    // ��� ����������� '.' ������ ','
    setlocale( LC_NUMERIC, "C" );


	//cout << endl << endl << endl << "----- testIntersectionCircleBox()" << endl << endl;
	//testIntersectionCircleBox();

	cout << endl << endl << endl << "----- testBurnTrianglePolygonAndTriangulateResult1()" << endl << endl;
	testBurnTrianglePolygonAndTriangulateResult1();

	//cout << endl << endl << endl << "----- testBurnTrianglePolygonAndTriangulateResult2()" << endl << endl;
	//testBurnTrianglePolygonAndTriangulateResult2();

	//cout << endl << endl << endl << "----- testBurnTrianglePolygonAndTriangulateResult3()" << endl << endl;
	//testBurnTrianglePolygonAndTriangulateResult3();


	cout << endl << "^" << endl;
    cin.ignore();

    return 0;
}
