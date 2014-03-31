#include "stdafx.h"
#include "external/triangle/triangle.h"


/**
* Радиус планеты, м.
*/
const double WORLD_PLANET_RADIUS = 5000.0 * 1000.0;

/**
* Сектор - квадратный участок, с которым работают визуализаторы.
* Задаётся в метрах.
*/
const double SECTOR_SIZE = 1000.0;
const double SECTOR_HALF_SIZE = SECTOR_SIZE / 2.0;


//typedef std::pair< double, double>  coord_xy_t;
typedef bgm::d2::point_xy< double >  p_t;
typedef bgm::polygon< p_t >          polygon_t;
typedef boost::tuple< p_t, p_t, p_t >  triangle_t;



/**
* @return Вершины круга. Представлены, по возможности, целыми числами.
*/
polygon_t drawCircle( double radius, size_t vertex ) {
	
	assert( (radius >= 100.00) && "Метод рассчитан на круг, радиусом >= 100.00" );

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
* Пересечение круга с прямоугольником. Круг намного больше прямоугольника.
*/
void testIntersectionCircleBox() {
    
    // Планета
	// (!) Слишком малая точность для круга может привести к тому, что
	//     сектор (см. ниже) и круг не пересекутся.
    const bgm::polygon< p_t >  circle = drawCircle( 5e6, 20 );
    cout << "Полигон планеты " << bg::dsv( circle ) << endl;
    cout << "Площадь " << bg::area( circle ) << endl << endl;


    // Сектор
    const auto aimCoord = p_t( 0.0, WORLD_PLANET_RADIUS );
    bgm::box< p_t >  sector(
        p_t( aimCoord.x() - SECTOR_HALF_SIZE,  aimCoord.y() - SECTOR_HALF_SIZE ),
        p_t( aimCoord.x() + SECTOR_HALF_SIZE,  aimCoord.y() + SECTOR_HALF_SIZE )
    );
    bg::correct( sector );
    cout << "Полигон сектора " << bg::dsv( sector ) << endl;
    cout << "Площадь " << bg::area( sector ) << endl << endl;


    // Пересечение
    std::vector< polygon_t > pv;
    bg::intersection( circle, sector, pv );
    cout << "Количество полигонов после пересечения " << pv.size() << endl << endl;
	assert( (pv.size() == 1) && "Ожидался 1 полигон." );

    // Анализируем полигон
    for (auto itr = pv.cbegin(); itr != pv.cend(); ++itr) {
        const polygon_t poly = *itr;
        // Каждый полигон - это отдельная реальная сущность
        cout << "Полигон планеты в секторе " << bg::dsv( poly ) << endl;
        cout << "Площадь " << bg::area( poly ) << endl << endl;
        // Вершины полигона
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
* Разбивает полигон на треугольники. Полигон может содержать дыры.
*
* @source http://www.cs.cmu.edu/~quake/triangle.html 
* @source external/triangle/tricall.c
* @source http://people.sc.fsu.edu/~jburkardt/c_src/triangle/triangle.html
* @see Описание параметров ниже, рядом с triangulate()
*
* Фигура описывается как набор всех точек. Дальше - два пути:
* 1. Фигура представляет собой цельное кольцо (без дыр).
*    В этом случае, формируем только список вершин 'pointlist'.
*
* 2. Фигура - кольцо с дырами.
*    Помимо 'pointlist' формируется:
*      а) список сегментов 'segmentlist'; представляет собой список вершин,
*         оформленных как кольца (всегда есть замыкающая вершина)
*         @example http://people.sc.fsu.edu/~jburkardt/data/triangle_files/double_hex2.poly
*                  (в примере внешнее кольцо задано 8-кой вершин для каждой стороны)
*      б) дыры; декларируются как координаты внутри заданных в 'segmentlist' колец
*
* @helper testBurnTrianglePolygonAndTriangulateResult()
*/
std::vector< triangle_t > triangulatePolygon( polygon_t& poly ) {

    // Полигон всегда завершаем координатами первой вершины
    bg::correct( poly );

    // Инициируем структуры, необходимые для работы метода триангуляции
    struct triangulateio in;
    struct triangulateio out;

    cout << "Полигон для триангуляции, внешнее кольцо: " <<
        poly.outer().size() << " вершины" << endl << endl;
    cout << "Полигон для триангуляции, внутреннИЕ кольцА: " <<
        poly.inners().size() << endl << endl;

	// Внешнее кольцо, всегда одно
    const auto outerVertex = poly.outer();
	// Последняя вершина дублирует первую, -1
	const auto countOuterVertex = outerVertex.size() - 1;

	// Внутренние кольца (дыры), могут отсутствовать или быть несколько
    const auto innerPoly = poly.inners();
	int countInnerVertex = 0;
	for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
        const auto hole = *itrPoly;
		// Последняя вершина дублирует первую, -1
		countInnerVertex += (int)hole.size() - 1;
	}


	const auto countVertex = countOuterVertex + countInnerVertex;

	cout << endl << "Заполняем 'pointlist'" << endl;

    in.numberofpoints = (int)countVertex;
    // Вершина задаётся парой координат
    in.pointlist = new real_t[ in.numberofpoints * 2 ];
	in.numberofpointattributes = 0;
    in.pointattributelist = nullptr;
    in.pointmarkerlist = nullptr;

    int pointerPoint2 = 0;

	// Внешнее кольцо (последняя вершина кольца дублирует первую, -1) ...
    for (auto itrVertex = outerVertex.cbegin(); itrVertex != outerVertex.cend() - 1; ++itrVertex) {
        const auto vertex = *itrVertex;
        const auto x = (real_t)vertex.x();
        const auto z = (real_t)vertex.y();
        cout << "Внешнее кольцо, вершина " << x << ", " << z << endl;
        in.pointlist[pointerPoint2 + 0] = x;
        in.pointlist[pointerPoint2 + 1] = z;
        pointerPoint2 += 2;
    }

	// ... плюс внутренние кольца
	for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
        const auto hole = *itrPoly;
		// (последняя вершина дыры дублирует первую, -1)
		for (auto itrVertex = hole.cbegin(); itrVertex != hole.cend() - 1; ++itrVertex) {
			const auto vertex = *itrVertex;
			const auto x = (real_t)vertex.x();
			const auto z = (real_t)vertex.y();
			cout << "Внутреннее кольцо, вершина " << x << ", " << z << endl;
			in.pointlist[pointerPoint2 + 0] = x;
			in.pointlist[pointerPoint2 + 1] = z;
			pointerPoint2 += 2;
		}

	} // for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly)


	// Сегменты и дыры заполняем только при наличия дыр
	if (countInnerVertex == 0) {

		in.numberofsegments = in.numberofholes = 0;
		in.segmentlist = in.segmentmarkerlist = nullptr;
		in.holelist = nullptr;

	} else {

		cout << endl << "Заполняем 'segmentlist'" << endl;

		// Кольца всегда замыкаются
		in.numberofsegments = (int)countVertex;
		in.segmentlist = new int[ in.numberofsegments * 2 ];
		in.segmentmarkerlist = nullptr;

		int pointerSegment2 = 0;
		int numberPoint = 0;
		int firstNumberPoint = numberPoint;

		// Внешнее кольцо и замыкающая вершина ...
		cout << endl << "Внешнее кольцо" << endl;
		for (auto itrVertex = outerVertex.cbegin(); itrVertex != outerVertex.cend() - 2; ++itrVertex) {
			in.segmentlist[pointerSegment2 + 0] = numberPoint;
			in.segmentlist[pointerSegment2 + 1] = numberPoint + 1;
			cout << "Внешнее кольцо, вершина " <<
				in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
			++numberPoint;
			pointerSegment2 += 2;
		}
		in.segmentlist[pointerSegment2 + 0] = numberPoint;
		in.segmentlist[pointerSegment2 + 1] = firstNumberPoint;
		cout << "Внешнее кольцо, замыкающая вершина " <<
			in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
		++numberPoint;
		pointerSegment2 += 2;

		// ... плюс внутренние кольца с замыкающими вершинами
		for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
			const auto hole = *itrPoly;
			cout << endl << "Внутреннее кольцо" << endl;
			firstNumberPoint = numberPoint;
			for (auto itrVertex = hole.cbegin(); itrVertex != hole.cend() - 2; ++itrVertex) {
				in.segmentlist[pointerSegment2 + 0] = numberPoint;
				in.segmentlist[pointerSegment2 + 1] = numberPoint + 1;
				cout << "Внутреннее кольцо, вершина " <<
					in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
				++numberPoint;
				pointerSegment2 += 2;
			}
			in.segmentlist[pointerSegment2 + 0] = numberPoint;
			in.segmentlist[pointerSegment2 + 1] = firstNumberPoint;
			cout << "Внутреннее кольцо, замыкающая вершина " <<
				in.segmentlist[pointerSegment2 + 0] << ", " << in.segmentlist[pointerSegment2 + 1] << endl;
			++numberPoint;
			pointerSegment2 += 2;

		} // for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly)


		// Указываем координатами, какие сегменты являются дырами (holes)
		in.numberofholes = (int)innerPoly.size();
		in.holelist = new real_t[ in.numberofholes * 2 ];

		int pointerHole2 = 0;

		for (auto itrPoly = innerPoly.cbegin(); itrPoly != innerPoly.cend(); ++itrPoly) {
			const auto hole = *itrPoly;
			p_t center;
			bg::centroid( hole, center );
			cout << endl << "Дыра, центр " << center.x() << ", " << center.y() << endl << endl;
			in.holelist[pointerHole2 + 0] = center.x();
			in.holelist[pointerHole2 + 1] = center.y();
			pointerHole2 += 2;
		}

	} // else if (countInnerVertex == 0)


	// Прочее
	in.numberofregions = 0;
	in.regionlist = nullptr;


    printf( "Input point set:\n\n" );
    printReportTriangulate( &in, 0, 0, 0, 1, 0, 0 );


    // Инициализируем структуры, в которые будет залит результат
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
	  * Звездой помечены полезные/интересные параметры.
	*/
	// @optimize '-X'
	// @optimize '-Q'
	// @optimize Что быстрее: '-i' или '-F'?
	// @interesting '-V'
    triangulate( "Xiz", &in, &out, nullptr );

    printf( "Initial triangulation:\n\n" );
    printReportTriangulate( &out, 0, 1, 0, 1, 0, 0 );

    // Собираем получившиеся треугольники
	// (!) triangle() не реагирует на параметр '-c' и всегда создаёт
	// из фигуры выпуклый многоугольник. Т.к. здесь нам нужен набор
	// треугольников, которые принадлежат исходной фигуре, проводим
	// в цикле их дополнительную проверку.
    assert( (out.trianglelist && (out.numberoftriangles > 0) ) && "Фигура не разбита на треугольники." );
    assert( (out.numberofcorners == 3) && "Ожидаем только треугольники." );
    std::vector< triangle_t >  vt;
    for (int i = 0; i < out.numberoftriangles; ++i) {
        const auto q = i * out.numberofcorners;

        // Получаем номер точки
        auto pn = out.trianglelist[ q + 0 ];
        // Получаем координаты точки
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

		// Проверяем, что треугольник принадлежит исходной фигуре
		// Проверку делаем через центр масс треугольника
		const auto c = bg::make< p_t >(
			(coord1.x() + coord2.x() + coord3.x()) / 3.0,
		    (coord1.y() + coord2.y() + coord3.y()) / 3.0
	    );
		if ( bg::within( c, poly ) ) {
			const triangle_t tri( coord1, coord2, coord3 );
			vt.push_back( tri );
		} else {
			cout << "Треугольник " <<
				bg::dsv( coord1 ) << " " <<
				bg::dsv( coord2 ) << " " <<
				bg::dsv( coord3 ) << " " << " исключён" <<
			endl << endl;
		}

    } // for (int i = 0; i < out.numberoftriangles; ++i)

    cout << "Собрано треугольников " << vt.size() << endl;


    // Убираем за собой
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
* Пересечение треугольника и кольца (полигон без дырок). Кольцо много
* меньше треугольника и лежит внутри треугольника, прожигая в нём дыру.
* После вычисления фигуры пересечения делаем её триангуляцию.
*/
void testBurnTrianglePolygonAndTriangulateResult1() {

	// Треугольник
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
	cout << "Треугольник: " << bg::dsv( tri ) << endl;

	// Кольцо внутри треугольника
	polygon_t ring;
	bg::append( ring, bg::make< p_t >( 10.0, 10.0 ) );
	bg::append( ring, bg::make< p_t >( 10.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 10.0 ) );
	bg::correct( ring );
	cout << "Кольцо внутри треугольника (дыра): " << bg::dsv( ring ) << endl;


	cout << endl << "Выжигаем в треугольнике дыру по форме кольца." << endl << endl;
	// @see http://www.boost.org/doc/libs/1_47_0/libs/geometry/doc/html/geometry/reference/algorithms/difference.html
    std::vector< polygon_t > figure;
    bg::difference( tri, ring, figure );
	
	cout << endl << "Полигонов после выжигания дыры: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "Полигон " << bg::dsv( *itr ) << endl;
    }

	assert( (figure.size() == 1) && "Должен получиться 1 полигон." );

	// Отдаём полученную фигуру на триангуляцию
	polygon_t firstPoly = *figure.cbegin();
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "Должно получиться много треугольников." );

}







/**
* @see testBurnTrianglePolygonAndTriangulateResult1()
*
* В треугольнике делается две неперекрывающиеся дыры.
*/
void testBurnTrianglePolygonAndTriangulateResult2() {

	// Треугольник
	polygon_t tri;
	bg::append( tri, bg::make< p_t >(   0.0,   0.0 ) );
	bg::append( tri, bg::make< p_t >(   0.0, 100.0 ) );
	bg::append( tri, bg::make< p_t >( 100.0,   0.0 ) );
	bg::correct( tri );
	cout << "Треугольник: " << bg::dsv( tri ) << endl;

	// Кольцо 1 внутри треугольника
	polygon_t ring1;
	bg::append( ring1, bg::make< p_t >( 10.0, 10.0 ) );
	bg::append( ring1, bg::make< p_t >( 10.0, 50.0 ) );
	bg::append( ring1, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring1, bg::make< p_t >( 20.0, 10.0 ) );
	bg::correct( ring1 );
	cout << "Кольцо 1 внутри треугольника (дыра): " << bg::dsv( ring1 ) << endl;

	// Кольцо 2 внутри треугольника, не перекрывает кольцо 1
	polygon_t ring2;
	bg::append( ring2, bg::make< p_t >( 1.0, 1.0 ) );
	bg::append( ring2, bg::make< p_t >( 1.0, 5.0 ) );
	bg::append( ring2, bg::make< p_t >( 2.0, 5.0 ) );
	bg::append( ring2, bg::make< p_t >( 2.0, 2.0 ) );
	bg::correct( ring2 );
	cout << "Кольцо 2 внутри треугольника (дыра): " << bg::dsv( ring2 ) << endl;
	assert( !bg::intersects( ring1, ring2) && "Кольца не должны пересекаться." );

	cout << endl << "Выжигаем в треугольнике дыру по форме кольца 1." << endl;
    std::vector< polygon_t > figure;
    bg::difference( tri, ring1, figure );
	assert( (figure.size() == 1) && "Должен получиться 1 полигон." );

	cout << endl << "Выжигаем в треугольнике дыру по форме кольца 2." << endl;
	polygon_t firstPoly = *figure.cbegin();
	figure.clear();
	bg::difference( firstPoly, ring2, figure );
	assert( (figure.size() == 1) && "Должен получиться 1 полигон." );
	
	cout << endl << "Полигонов после выжигания 2 дыр: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "Полигон " << bg::dsv( *itr ) << endl;
    }

	// Отдаём полученную фигуру на триангуляцию
	firstPoly = *figure.cbegin();
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "Должно получиться много треугольников." );

}







/**
* @see testBurnTrianglePolygonAndTriangulateResult1()
*
* В треугольнике делается дыра, делающая из него вогнутый полигон.
*/
void testBurnTrianglePolygonAndTriangulateResult3() {

	// Треугольник
	polygon_t tri;
	bg::append( tri, bg::make< p_t >(   0.0,   0.0 ) );
	bg::append( tri, bg::make< p_t >(   0.0, 100.0 ) );
	bg::append( tri, bg::make< p_t >( 100.0,   0.0 ) );
	bg::correct( tri );
	cout << "Треугольник: " << bg::dsv( tri ) << endl;

	// Кольцо, часть - внутри треугольника, часть - снаружи
	polygon_t ring;
	bg::append( ring, bg::make< p_t >( -5.0, -5.0 ) );
	bg::append( ring, bg::make< p_t >( -5.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, 50.0 ) );
	bg::append( ring, bg::make< p_t >( 20.0, -5.0 ) );
	bg::correct( ring );
	cout << "Кольцо выбивает часть треугольника: " << bg::dsv( ring ) << endl;
	//assert( bg::overlaps( tri, ring ) && "Треугольник и кольцо должны частично перекрывать друг друга." );
	//assert( bg::intersects( tri, ring ) && !bg::within( ring, tri ) && "Треугольник и кольцо должны частично перекрывать друг друга." );


	cout << endl << "Выбиваем в треугольнике часть по форме кольца." << endl << endl;
    std::vector< polygon_t > figure;
    bg::difference( tri, ring, figure );
	
	cout << endl << "Полигонов после выбивания: " << figure.size() << endl;
	for (auto itr = figure.cbegin(); itr != figure.cend(); ++itr) {
	    cout << "Полигон " << bg::dsv( *itr ) << endl;
    }

	assert( (figure.size() == 1) && "Должен получиться 1 полигон." );

	// Отдаём вогнутый полигон на триангуляцию
	polygon_t firstPoly = *figure.cbegin();
	assert( (firstPoly.inners().size() == 0) && "Фигура должна быть с выбоиной, но без дыры." );
	const std::vector< triangle_t >  r = triangulatePolygon( firstPoly );
	assert( (r.size() > 1) && "Должно получиться много треугольников." );

}







int main() {

	setlocale( LC_ALL, "Russian" );
    // Для разделителя '.' вместо ','
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
