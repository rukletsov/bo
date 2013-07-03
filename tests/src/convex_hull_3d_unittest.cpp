
#include <gtest/gtest.h>

#include "bo/surfaces/convex_hull_3d.hpp"

using namespace bo::surfaces;

// "Test fixture" using base class form GTEST.
class ConvexHull3DTest: public testing::Test
{
protected:
    typedef IncrementalConvexHull3D<float> Hull;
    typedef Hull::Point3D Point3D;
    typedef Hull::Points3D Points3D;

    ConvexHull3DTest()
    { }

    virtual void SetUp()
    {
        p_small_.push_back(Point3D(0, 0, 0));
        p_small_.push_back(Point3D(2.71f, -3.14f, 5.27f));
        p_small_.push_back(Point3D(100, -1.5f, -8));

        p_collinear_.push_back(Point3D(0, 0, 0));
        p_collinear_.push_back(Point3D(0.5f, 0.5f, 0.5f));
        p_collinear_.push_back(Point3D(1.2f, 1.2f, 1.2f));
        p_collinear_.push_back(Point3D(3.14f, 3.14f, 3.14f));

        p_coplanar_.push_back(Point3D(0, 0, 0));
        p_coplanar_.push_back(Point3D(1, 0, 0));
        p_coplanar_.push_back(Point3D(0, 1, 0));
        p_coplanar_.push_back(Point3D(0.5f, 0.5f, 0.00001f));

        p_unit_cube_.push_back(Point3D(0, 0, 0));
        p_unit_cube_.push_back(Point3D(1, 0, 0));
        p_unit_cube_.push_back(Point3D(0, 1, 0));
        p_unit_cube_.push_back(Point3D(0, 0, 1));
        p_unit_cube_.push_back(Point3D(1, 1, 0));
        p_unit_cube_.push_back(Point3D(0, 1, 1));
        p_unit_cube_.push_back(Point3D(1, 0, 1));
        p_unit_cube_.push_back(Point3D(1, 1, 1));
        p_unit_cube_.push_back(Point3D(0.5f, 0.5f, 0.5f));

        p_central_cube_.push_back(Point3D(-0.25f, -0.25f, -0.25f));
        p_central_cube_.push_back(Point3D(0.25f, -0.25f, -0.25f));
        p_central_cube_.push_back(Point3D(-0.25f, 0.25f, -0.25f));
        p_central_cube_.push_back(Point3D(-0.25f, -0.25f, 0.25f));
        p_central_cube_.push_back(Point3D(0.25f, 0.25f, -0.25f));
        p_central_cube_.push_back(Point3D(-0.25f, 0.25f, 0.25f));
        p_central_cube_.push_back(Point3D(0.25f, -0.25f, 0.25f));
        p_central_cube_.push_back(Point3D(0.25f, 0.25f, 0.25f));
        p_central_cube_.push_back(Point3D(0.2f, -0.2f, 0.2f));

        p_small_faces_.push_back(Point3D(1.20265543f, 2.7372396f, -1.91974759f));
        p_small_faces_.push_back(Point3D(1.46370459f, 1.4345957f, -1.34324634f));
        p_small_faces_.push_back(Point3D(1.82530951f, 2.51767373f, -3.18732214f));
        p_small_faces_.push_back(Point3D(2.44798946f, 3.33310866f, -1.60480368f));
        p_small_faces_.push_back(Point3D(1.8206923f, 2.52521658f, -3.1833353f));
        p_small_faces_.push_back(Point3D(2.08507633f, 1.22142851f, -2.61365271f));
        p_small_faces_.push_back(Point3D(2.70903873f, 2.03046465f, -1.02830243f));
        p_small_faces_.push_back(Point3D(1.82530951f, 2.51767373f, -3.18732214f));
        p_small_faces_.push_back(Point3D(3.29080701f, 2.01492167f, -2.38616967f));
        p_small_faces_.push_back(Point3D(3.33041072f, 1.8172977f, -2.29870892f));

        p_parallel_face_.push_back(Point3D(2.93130445f, 0.924954236f, -1.66728044f));
        p_parallel_face_.push_back(Point3D(3.23396015f, 0.36961019f, -1.74894285f));
        p_parallel_face_.push_back(Point3D(3.4822154f, 0.828871608f, -2.43442965f));
        p_parallel_face_.push_back(Point3D(4.45183563f, 0.477457225f, -0.887860775f));
        p_parallel_face_.push_back(Point3D(3.42953944f, 1.84666657f, -0.869312644f));
        p_parallel_face_.push_back(Point3D(3.79619265f, 2.06103563f, -0.968243122f));
        p_parallel_face_.push_back(Point3D(4.09512568f, 1.18721819f, -2.59980488f));
        p_parallel_face_.push_back(Point3D(3.69292235f, 1.86999011f, -0.683091283f));
    }

    Points3D p_small_;
    Points3D p_collinear_;
    Points3D p_coplanar_;
    Points3D p_unit_cube_;
    Points3D p_central_cube_;
    Points3D p_small_faces_;
    Points3D p_parallel_face_;
};

TEST_F(ConvexHull3DTest, SmallCloud)
{
    Hull conhul3d(p_small_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 0U);

    float volume = conhul3d.get_volume();

    EXPECT_EQ(volume, 0);
}

TEST_F(ConvexHull3DTest, Collinear)
{
    Hull conhul3d(p_collinear_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 0U);

    float volume = conhul3d.get_volume();

    EXPECT_EQ(volume, 0);
}

TEST_F(ConvexHull3DTest, Coplanar)
{
    Hull conhul3d(p_coplanar_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 0U);

    float volume = conhul3d.get_volume();

    EXPECT_EQ(volume, 0);
}

TEST_F(ConvexHull3DTest, UnitCube)
{
    Hull conhul3d(p_unit_cube_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 12U);

    float volume = conhul3d.get_volume();

    EXPECT_NEAR(volume, 1.0f, 0.001f);
}

TEST_F(ConvexHull3DTest, CentralCube)
{
    Hull conhul3d(p_central_cube_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 12U);

    float volume = conhul3d.get_volume();

    EXPECT_NEAR(volume, 0.125f, 0.001f);
}

TEST_F(ConvexHull3DTest, SmallFacesVolume)
{
    Hull conhul3d(p_small_faces_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 14U);

    float volume = conhul3d.get_volume();

    EXPECT_NEAR(volume, 2.45f, 0.01f);
}

TEST_F(ConvexHull3DTest, AddingParallelFace)
{
    Hull conhul3d(p_parallel_face_);
    Hull::Faces3D f = conhul3d.get_faces();

    EXPECT_EQ(f.size(), 12U);

    float volume = conhul3d.get_volume();

    EXPECT_NEAR(volume, 1.0f, 0.01f);
}
