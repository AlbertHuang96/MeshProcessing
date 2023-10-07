// MeshDelaunay.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <openmesh/Core/IO/MeshIO.hh>
#include "Eigen/Dense"

#include <queue>
#include <iostream>

typedef OpenMesh::TriMesh_ArrayKernelT<> MeshType;
OpenMesh::TriMesh_ArrayKernelT<> mesh2;

auto Vertices2 = OpenMesh::makeTemporaryProperty<MeshType::VertexHandle, Eigen::Vector3f>(mesh2);

//--------------------------Delaunay triangulation--------------------------//

double getTriFaceArea(MeshType::VertexFaceIter vfIter)
{
	MeshType::VertexHandle vh0 = vfIter->halfedge().from();
	MeshType::VertexHandle vh1 = vfIter->halfedge().to();
	MeshType::VertexHandle vh2 = vfIter->halfedge().next().to();

	Eigen::Vector3f v0 = Vertices2[vh0];
	Eigen::Vector3f v1 = Vertices2[vh1];
	Eigen::Vector3f v2 = Vertices2[vh2];

	/*Eigen::Vector3f v0, v1, v2;
	v0[0] = mesh2.point(vh0)[0];
	v0[1] = mesh2.point(vh0)[1];
	v0[2] = mesh2.point(vh0)[2];

	v1[0] = mesh2.point(vh1)[0];
	v1[1] = mesh2.point(vh1)[1];
	v1[2] = mesh2.point(vh1)[2];

	v2[0] = mesh2.point(vh2)[0];
	v2[1] = mesh2.point(vh2)[1];
	v2[2] = mesh2.point(vh2)[2];*/

	Eigen::Vector3f a = v1 - v0;
	Eigen::Vector3f b = v2 - v0;
	Eigen::Vector3f area(a[1] * b[2] - b[1] * a[2], -(a[0] * b[2] - b[0] * a[2]), a[0] * b[1] - b[0] * a[1]);
	return 0.5 * area.norm();
}

Eigen::Vector3f getTriFaceCircumcenter(MeshType::VertexFaceIter vfIter)
{
	MeshType::VertexHandle vh0 = vfIter->halfedge().from();
	MeshType::VertexHandle vh1 = vfIter->halfedge().to();
	MeshType::VertexHandle vh2 = vfIter->halfedge().next().to();
	Eigen::Vector3f v0 = Vertices2[vh0];
	Eigen::Vector3f v1 = Vertices2[vh1];
	Eigen::Vector3f v2 = Vertices2[vh2];

	/*Eigen::Vector3f v0, v1, v2;
	v0[0] = mesh2.point(vh0)[0];
	v0[1] = mesh2.point(vh0)[1];
	v0[2] = mesh2.point(vh0)[2];

	v1[0] = mesh2.point(vh1)[0];
	v1[1] = mesh2.point(vh1)[1];
	v1[2] = mesh2.point(vh1)[2];

	v2[0] = mesh2.point(vh2)[0];
	v2[1] = mesh2.point(vh2)[1];
	v2[2] = mesh2.point(vh2)[2];*/

	double x1, y1, x2, y2, x3, y3;
	x1 = v0[0]; y1 = v0[1];
	x2 = v1[0]; y2 = v1[1];
	x3 = v2[0]; y3 = v2[1];

	double a1, b1, c1, a2, b2, c2;
	a1 = 2 * (x2 - x1);	  a2 = 2 * (x3 - x2);	c1 = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
	b1 = 2 * (y2 - y1);	  b2 = 2 * (y3 - y2);	c2 = x3 * x3 + y3 * y3 - x2 * x2 - y2 * y2;

	Eigen::Vector3f circumcenter(0.0, 0.0, 0.0);
	circumcenter[0] = (b2 * c1 - b1 * c2) / (a1 * b2 - a2 * b1);
	circumcenter[1] = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);
	circumcenter[2] = 0;

	return circumcenter;
}

void DelaunayTriangulation(int numIter)
{
	//OpenMesh::IO::read_mesh(mesh, "models/venus.obj");
	//OpenMesh::IO::read_mesh(mesh, "models/chinesedragon.obj");
	if (!OpenMesh::IO::read_mesh(mesh2, "models/leaf.obj"))
	{
		return;
	}

	clock_t start, end;
	start = clock();

	if (mesh2.n_vertices() == 3)
	{
		return;
	}
	int numVertices = mesh2.n_vertices();

	for (MeshType::VertexIter vIter = mesh2.vertices_begin(); vIter != mesh2.vertices_end(); ++vIter)
	{
		Vertices2[*vIter][0] = mesh2.point(*vIter)[0];
		Vertices2[*vIter][1] = mesh2.point(*vIter)[1];
		Vertices2[*vIter][2] = mesh2.point(*vIter)[2];
	}

	for (int i = 0; i < numIter; i++)
	{
		//flip
		for (MeshType::EdgeIter eIter = mesh2.edges_begin(); eIter != mesh2.edges_end(); ++eIter)
		{
			if (mesh2.is_boundary(*eIter) || !mesh2.is_flip_ok(*eIter))
			{
				continue;
			}

			MeshType::VertexHandle vh1 = eIter->halfedge(0).from();
			MeshType::VertexHandle vh2 = eIter->halfedge(0).to();
			MeshType::VertexHandle vh3 = eIter->halfedge(0).next().to();
			MeshType::VertexHandle vh4 = eIter->halfedge(1).next().to();
			Eigen::Vector3f v1 = Vertices2[vh1];
			Eigen::Vector3f v2 = Vertices2[vh2];
			Eigen::Vector3f v3 = Vertices2[vh3];
			Eigen::Vector3f v4 = Vertices2[vh4];

			/*Eigen::Vector3f v1, v2, v3, v4;

			v1[0] = mesh2.point(vh1)[0];
			v1[1] = mesh2.point(vh1)[1];
			v1[2] = mesh2.point(vh1)[2];

			v2[0] = mesh2.point(vh2)[0];
			v2[1] = mesh2.point(vh2)[1];
			v2[2] = mesh2.point(vh2)[2];

			v3[0] = mesh2.point(vh3)[0];
			v3[1] = mesh2.point(vh3)[1];
			v3[2] = mesh2.point(vh3)[2];

			v4[0] = mesh2.point(vh4)[0];
			v4[1] = mesh2.point(vh4)[1];
			v4[2] = mesh2.point(vh4)[2];*/

			double alpha(0.0), alpha1(0.0), alpha2(0.0);
			alpha1 = acos((pow((v1 - v3).norm(), 2) + pow((v2 - v3).norm(), 2)
				- pow((v1 - v2).norm(), 2)) / (2 * (v1 - v3).norm() * (v2 - v3).norm()));
			alpha2 = acos((pow((v1 - v4).norm(), 2) + pow((v2 - v4).norm(), 2)
				- pow((v1 - v2).norm(), 2)) / (2 * (v1 - v4).norm() * (v2 - v4).norm()));
			alpha = alpha1 + alpha2;
			if (alpha > M_PI)
			{
				mesh2.flip(*eIter);
			}
		}

		//update vertices
		for (MeshType::VertexIter vIter = mesh2.vertices_begin(); vIter != mesh2.vertices_end(); ++vIter)
		{
			if (mesh2.is_boundary(*vIter))
			{
				continue;
			}
			Eigen::Vector3f tmp(0.0f, 0.0f, 0.0f);
			double area(0.0), sum_area(0.0);
			for (MeshType::VertexFaceIter vfIter = mesh2.vf_iter(*vIter); vfIter.is_valid(); ++vfIter)
			{
				area = getTriFaceArea(vfIter);
				sum_area += area;
				Eigen::Vector3f center = getTriFaceCircumcenter(vfIter);
				tmp = tmp + area * center;
			}

			mesh2.set_point(*vIter, MeshType::Point(tmp[0] / sum_area, tmp[1] / sum_area, tmp[2] / sum_area));
		}
	}

	end = clock();
	std::cout << "time needed to triangulate = " << (double)(end - start) / CLOCKS_PER_SEC << "second" << std::endl;
	// write mesh
	OpenMesh::IO::write_mesh(mesh2, "models/leaf3.obj");
}

int main()
{
	DelaunayTriangulation(1000);
    //std::cout << "Hello World!\n";
}

