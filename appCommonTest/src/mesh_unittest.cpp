
#include <cmath>
#include <gtest/gtest.h>

#include "common/mesh.hpp"
#include "debug_alloc.hpp"

using namespace common;

// Directory where test data is stored. This variable is set during the configuration
// of the test environment and should be used to determine where external data for
// the tests is stored.
extern std::string DataDirectory;


// Create a so-called "text fixture" using base class form GTEST.
class MeshTest: public testing::Test
{
protected:

    MeshTest(): mesh1_(5)
    { }

    virtual void SetUp()
    {
        // Add 5 vertices and remember their indices.
        vertex1_ = mesh1_.add_vertex(Mesh::Vertex(1., 0., 1.));
        vertex2_ = mesh1_.add_vertex(Mesh::Vertex(2., 1., 1.));
        vertex3_ = mesh1_.add_vertex(Mesh::Vertex(3., 0., 1.));
        vertex4_ = mesh1_.add_vertex(Mesh::Vertex(2., -1., 1.));
        vertex5_ = mesh1_.add_vertex(Mesh::Vertex(2., 0., 2.));

        // Specify faces in the mesh using previously remembered vertices' indices.
        basis1_ = mesh1_.add_face(Mesh::Face(vertex1_, vertex2_, vertex4_));
        basis2_ = mesh1_.add_face(Mesh::Face(vertex2_, vertex3_, vertex4_));
        slope1_ = mesh1_.add_face(Mesh::Face(vertex1_, vertex2_, vertex5_));
        mesh1_.add_face(Mesh::Face(vertex1_, vertex4_, vertex5_));
        mesh1_.add_face(Mesh::Face(vertex3_, vertex2_, vertex5_));
        mesh1_.add_face(Mesh::Face(vertex3_, vertex4_, vertex5_));
    }

    // Create a test mesh object with 5 vertices: a square pyramid with the edge
    // lenght equals sqrt(2). Its basis is parallel to OXY surface (z = 1).
    Mesh mesh1_;

    // Indices of the certain faces and vertices.
    std::size_t vertex1_;
    std::size_t vertex2_;
    std::size_t vertex3_;
    std::size_t vertex4_;
    std::size_t vertex5_;
    std::size_t basis1_;
    std::size_t basis2_;
    std::size_t slope1_;
};


TEST_F(MeshTest, SetUp)
{
    // Check if the pyramid mesh object was created successfully.
    ASSERT_EQ(std::size_t(5), mesh1_.get_all_vertices().size());
    ASSERT_EQ(std::size_t(6), mesh1_.get_all_faces().size());
    ASSERT_EQ(std::size_t(6), mesh1_.get_all_face_normals().size());

    // The normal direction can differ, but the normals should be of the same
    // length and collinear.
    ASSERT_DOUBLE_EQ(1., abs(
        (Mesh::Normal(0., 0., 1.) * mesh1_.get_all_face_normals()[basis1_])));

    ASSERT_DOUBLE_EQ(1., abs(
        (Mesh::Normal(0., 0., 1.) * mesh1_.get_all_face_normals()[basis2_])));

    ASSERT_DOUBLE_EQ(sqrt(3.) / 3., abs(
        (Mesh::Normal(0., 0., 1.) * mesh1_.get_all_face_normals()[slope1_])));
}
