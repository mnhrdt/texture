#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <fstream>
#include <map>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <vector>

#include "pickopt.c"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef Polyhedron::Vertex_handle   Vertex_handle;

typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};

int main_remesh(int argc, char* argv[])
{
    if (argc < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s input_mesh.off output_mesh.off\n",*argv);
  const char* filename_in = argv[1];
  const char* filename_out = argv[2];
  std::ifstream input(filename_in);
  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }
  double target_edge_length = 0.45;
  unsigned int nb_iter = 10;
  std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh),
      mesh,
      boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
  std::cout << "done." << std::endl;
  std::cout << "Start remeshing of " << filename_in
    << " (" << num_faces(mesh) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      .protect_constraints(true)//i.e. protect border, here
      );
  std::cout << "Remeshing done." << std::endl;
  std::ofstream remeshed(filename_out);
  remeshed << mesh;
  remeshed.close();
  return 0;
}

int main_refine(int c, char* v[])
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
            CGAL::Polygon_mesh_processing::parameters::density_control_factor(22.));
    //std::ofstream refined_off("refined.off");
    std::ofstream refined_off(filename_out);
    refined_off << poly;
    refined_off.close();
    std::cout << "Refinement added " << new_vertices.size() << " vertices." << std::endl;
    return 0;
}

int main(int c, char *v[])
{
    int method = atoi(pick_option(&c, &v, "m", "0"));
    if (method == 0)
        return main_remesh(c,v);
    else
        return main_refine(c,v);
}
