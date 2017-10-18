#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <fstream>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef Polyhedron::Vertex_handle   Vertex_handle;


int main(int c, char* v[])
{
    if (c < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s input_mesh.off output_mesh.off\n",*v);
    const char* filename_in = v[1];
    const char* filename_out = v[2];
    std::ifstream input(filename_in);
    Polyhedron poly;
    if ( !input || !(input >> poly) || poly.empty()
            || !CGAL::is_triangle_mesh(poly)) {
        std::cerr << "Not a valid input file." << std::endl;
        return 1;
    }
    std::vector<Polyhedron::Facet_handle>  new_facets;
    std::vector<Vertex_handle> new_vertices;
    CGAL::Polygon_mesh_processing::refine(poly,
            faces(poly),
            std::back_inserter(new_facets),
            std::back_inserter(new_vertices),
            CGAL::Polygon_mesh_processing::parameters::density_control_factor(2.));
    //std::ofstream refined_off("refined.off");
    std::ofstream refined_off(filename_out);
    refined_off << poly;
    refined_off.close();
    std::cout << "Refinement added " << new_vertices.size() << " vertices." << std::endl;
    return 0;
}

