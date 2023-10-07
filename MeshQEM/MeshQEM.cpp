// MeshQEM.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <queue>

#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <openmesh/Core/IO/MeshIO.hh>
#include "Eigen/Dense"


typedef OpenMesh::TriMesh_ArrayKernelT<> MeshType;
OpenMesh::TriMesh_ArrayKernelT<> mesh;

auto Qv = OpenMesh::makeTemporaryProperty<MeshType::VertexHandle, Eigen::Matrix4f>(mesh);
auto State = OpenMesh::makeTemporaryProperty<MeshType::EdgeHandle, int>(mesh);
auto Vertices = OpenMesh::makeTemporaryProperty<MeshType::VertexHandle, Eigen::Vector4f>(mesh);

//Eigen::Vector4f [a, b, c, d] for ax + by + cz + d = 0
auto planeParams = OpenMesh::makeTemporaryProperty<MeshType::FaceHandle, Eigen::Vector4f>(mesh);

struct edgeCollapse
{
	int state;
	float cost;
	MeshType::Point point;
	MeshType::HalfedgeHandle halfEdge;
	MeshType::EdgeHandle edge;
	MeshType::VertexHandle vertexTo;
	MeshType::VertexHandle vertexFrom;
};

struct cmp
{
	bool operator()(const edgeCollapse& a, const edgeCollapse& b)
	{
		return a.cost > b.cost;
	}
};

std::priority_queue<edgeCollapse, std::vector<edgeCollapse>, cmp> Cost;

//--------------------------mesh simplify--------------------------//
void computeQForVertex(MeshType::VertexHandle vh, MeshType& mesh)
{
	if (!mesh.is_boundary(vh))
	{
		Eigen::Matrix4f mat;

		for (MeshType::VertexFaceIter vfIter = mesh.vf_iter(vh); vfIter.is_valid(); ++vfIter)
		{
			mat += planeParams[*vfIter] * planeParams[*vfIter].transpose();
		}
		Qv[vh] = mat;
	}
}

void computeCost(MeshType::EdgeIter eIter)
{
	Eigen::Matrix4f tempQ = Qv[eIter->v0()] + Qv[eIter->v1()];

	Eigen::Matrix4f Q = tempQ;
	Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
	Q(3, 0) = 0.0f;
	Q(3, 1) = 0.0f;
	Q(3, 2) = 0.0f;
	Q(3, 3) = 1.0f;
	Eigen::FullPivLU<Eigen::Matrix4f> lu(Q);
	Eigen::Vector4f newPoint;
	if (lu.isInvertible())
	{
		newPoint = Q.inverse() * b;
	}
	else
	{
		newPoint = (Vertices[eIter->v0()] + Vertices[eIter->v1()]) / 2;
	}
	edgeCollapse eCollapse;
	eCollapse.edge = *eIter;
	eCollapse.halfEdge = eIter->halfedge(0);
	eCollapse.vertexTo = eIter->halfedge(0).to();
	eCollapse.vertexFrom = eIter->halfedge(0).from();
	eCollapse.cost = newPoint.transpose() * tempQ * newPoint;
	eCollapse.point = MeshType::Point(newPoint[0], newPoint[1], newPoint[2]);
	eCollapse.state = State[*eIter];

	Cost.push(eCollapse);
}

void updateCost(MeshType::VertexHandle vh, MeshType& mesh)
{
	for (MeshType::VertexOHalfedgeIter hfIter = mesh.voh_iter(vh); hfIter.is_valid(); ++hfIter)
	{
		MeshType::EdgeHandle edge = hfIter->edge();
		State[edge]++;
		MeshType::VertexHandle toVertex = hfIter->to();
		Eigen::Matrix4f newQ = Qv[vh] + Qv[toVertex];
		Eigen::Matrix4f tempQ = newQ;
		Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
		tempQ(3, 0) = 0.0f;
		tempQ(3, 1) = 0.0f;
		tempQ(3, 2) = 0.0f;
		tempQ(3, 3) = 1.0f;
		Eigen::FullPivLU<Eigen::Matrix4f> lu(tempQ);
		Eigen::Vector4f newPoint;
		if (lu.isInvertible())
		{
			newPoint = tempQ.inverse() * b;
		}
		else
		{
			newPoint = (Vertices[vh] + Vertices[toVertex]) / 2.0f;
		}

		edgeCollapse eCollapse;
		eCollapse.edge = edge;
		eCollapse.halfEdge = *hfIter;
		eCollapse.cost = newPoint.transpose() * newQ * newPoint;
		eCollapse.point = MeshType::Point(newPoint[0], newPoint[1], newPoint[2]);
		eCollapse.vertexFrom = vh;
		eCollapse.vertexTo = toVertex;
		eCollapse.state = State[edge];
		Cost.push(eCollapse);
	}
}

void updateQv(MeshType::VertexHandle vh, MeshType& mesh)
{
	computeQForVertex(vh, mesh);

	//update 1-ring neighboring vertices
	for (auto vvIter = mesh.vv_begin(vh); vvIter.is_valid(); ++vvIter)
	{
		computeQForVertex(*vvIter, mesh);
	}
}

bool collapse(const edgeCollapse& eCollapse, MeshType& mesh)
{
	bool bCollapse = false;
	MeshType::VertexHandle vh;
	if (mesh.is_collapse_ok(eCollapse.halfEdge))
	{
		mesh.collapse(eCollapse.halfEdge);
		vh = eCollapse.vertexTo;
		bCollapse = true;
	}
	else if (mesh.is_collapse_ok(mesh.opposite_halfedge_handle(eCollapse.halfEdge)))
	{
		mesh.collapse(mesh.opposite_halfedge_handle(eCollapse.halfEdge));
		vh = eCollapse.vertexFrom;
		bCollapse = true;
	}
	if (bCollapse)
	{
		updateQv(vh, mesh);

		mesh.set_point(vh, eCollapse.point);
		Vertices[vh][0] = eCollapse.point[0];
		Vertices[vh][1] = eCollapse.point[1];
		Vertices[vh][2] = eCollapse.point[2];
		Vertices[vh][0] = 1.0f;

		updateCost(vh, mesh);
	}

	return bCollapse;
}

void simplifyMesh(float ratio)
{
	//OpenMesh::IO::read_mesh(mesh, "models/venus.obj");
	//OpenMesh::IO::read_mesh(mesh, "models/chinesedragon.obj");
	if (!OpenMesh::IO::read_mesh(mesh, "models/leaf.obj"))
	{
		return;
	}

	std::cout << "number of the faces of this mesh = " << mesh.n_faces() << std::endl;

	clock_t start, end;
	start = clock();

	if (mesh.n_vertices() == 3)
	{
		return;
	}

	int numVertices = mesh.n_vertices();
	int numTargets = (1.0f - ratio) * mesh.n_vertices();

	mesh.request_face_normals();

	// calculate planes' parameters for Qv
	for (MeshType::FaceIter ft = mesh.faces_begin(); ft != mesh.faces_end(); ++ft)
	{
		float a, b, c, d;

		a = mesh.normal(*ft)[0];
		b = mesh.normal(*ft)[1];
		c = mesh.normal(*ft)[2];
		MeshType::Point pt = mesh.point(ft->halfedge().to());
		d = -(a * pt[0] + b * pt[1] + c * pt[2]);
		planeParams[*ft][0] = a;
		planeParams[*ft][1] = b;
		planeParams[*ft][2] = c;
		planeParams[*ft][3] = d;
	}

	//initial Qv
	for (MeshType::VertexIter vIter = mesh.vertices_begin(); vIter != mesh.vertices_end(); ++vIter)
	{
		computeQForVertex(*vIter, mesh);
		Vertices[*vIter][0] = mesh.point(*vIter)[0];
		Vertices[*vIter][1] = mesh.point(*vIter)[1];
		Vertices[*vIter][2] = mesh.point(*vIter)[2];
		Vertices[*vIter][3] = 1.0f;
	}

	//initial cost
	for (MeshType::EdgeIter eIter = mesh.edges_begin(); eIter != mesh.edges_end(); ++eIter)
	{
		State[*eIter] = 0;
		computeCost(eIter);
	}

	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	// start to collapse
	int i = 0;
	while (i < numTargets)
	{
		edgeCollapse eCollapse = Cost.top();
		Cost.pop();
		if (mesh.status(eCollapse.vertexFrom).deleted() || mesh.status(eCollapse.vertexTo).deleted())
		{
			continue;
		}

		/*MeshType::EdgeHandle edge;
		for (MeshType::EdgeIter eIter = mesh.edges_begin(); eIter != mesh.edges_end(); ++eIter)
		{
			if (eIter->halfedge(0) == eCollapse.halfEdge || eIter->halfedge(1) == eCollapse.halfEdge)
			{
				edge = eIter;
			}
		}*/

		if (eCollapse.state == State[eCollapse.edge])
		{
			if (collapse(eCollapse, mesh))
			{
				i++;
			}
		}
	}

	end = clock();
	std::cout << "time needed to collapse = " << (double)(end - start) / CLOCKS_PER_SEC << "second" << std::endl;
	// write mesh
	OpenMesh::IO::write_mesh(mesh, "models/leaf4.obj");
}

int main()
{
	simplifyMesh(0.5);
    //std::cout << "Hello World!\n";
}

