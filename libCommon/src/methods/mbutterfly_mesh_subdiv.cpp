#include "common/methods/mbutterfly_mesh_subdiv.hpp"
#include <map>
#include <list>

namespace common
{

namespace methods
{

Mesh::Vertex getAveragedVertex(const Mesh &m, size_t inda, size_t indb)
{
	Mesh::Faces f = m.get_all_faces();
	Mesh::Vertices v = m.get_all_vertices();

	Mesh::AdjacentFacesPerVertex adjacentFaces = m.get_neighbouring_faces_by_vertex(inda);
	Mesh::AdjacentVerticesPerVertex adjacentVertices = m.get_neighbouring_vertices(inda);

	//If the central vertex is on the boundary, calculate
	//the average value using the boundary schema
	if(adjacentFaces.size() != adjacentVertices.size())
	{
		//return v[inda]/2;
		return Mesh::Vertex(0,0,0);
	}
		
	//Array of the sorted adjacent vertices (clockwise or counter clockwise)
	std::vector<size_t> svert;

	//Sort the adjacent vertices(edges)
	{	
		size_t cur = indb;

		while(adjacentFaces.size()>0)
		{
			//Walk through the set of the adjacent faces of inda consequently in 
			//(counter) clockwise order beginning from the face that contains intb,
			//and save the resulting order of the vertices.
			Mesh::AdjacentFacesPerVertex::iterator it = adjacentFaces.begin();
			while(it != adjacentFaces.end())
			{
				if(f[*it].A() == cur || f[*it].B() == cur || f[*it].C() == cur)
				{
					//Save the order
					svert.push_back(cur);

					//Find among three vertices of the given triangle face
					//the vertex that is not equal neither inda nor cur, 
					//and update the value of cur with it.
					if((f[*it].A() == inda && f[*it].B() == cur) || (f[*it].B() == inda && f[*it].A() == cur))
						cur=f[*it].C();
					else 
					if((f[*it].B() == inda && f[*it].C() == cur) || (f[*it].C() == inda && f[*it].B() == cur))
						cur=f[*it].A();
					else 
					if((f[*it].C() == inda && f[*it].A() == cur) || (f[*it].A() == inda && f[*it].C() == cur))
						cur=f[*it].B();

					//Remove the processed face 
					adjacentFaces.erase(it);
					break;
				}

				++it;
			}
		}
	}

	//Get the valence of the stencil center
	int k = svert.size();

	//Realize the Modified Butterfly subdivision surface scheme
	Mesh::Vertex average(0,0,0);
	{   
        //Regular case
        if(k == 6) 
		{
			float a, b, c, d;
			a = 1/2.0f;
			b = 1/8.0f;
			c = -1/16.0f;
			d = 0.0f;
			
			average = a*v[svert[0]] + b/2*v[svert[1]] + c*v[svert[2]]
					  + d*v[svert[3]] + c*v[svert[4]] + b/2*v[svert[5]];
		}
		//Non-regular cases
		else if(k >= 5)
		{
			float pi = 3.14159f;
			for(int j = 0; j < k; ++j)
			{
				float s_j = (0.25f + cos(2*pi*j/k) + 0.5f*cos(4*pi*j/k))/k;
				average += s_j*v[svert[j]];
			}
		}
		else if(k == 3)
		{
			float s_0, s_1, s_2;
			s_0 = 5/12.0f;
			s_1 = -1/12.0f;
			s_2 = s_1;
			average = s_0*v[svert[0]] + s_1*v[svert[1]] + s_2*v[svert[2]];
		}
		else if(k == 4)
		{
			float s_0, s_1, s_2, s_3;
			s_0 = 3/8.0f;
			s_1 = 0;
			s_2 = -1/8.0f;
			s_3 = 0;
			average = s_0*v[svert[0]] + s_1*v[svert[1]] + s_2*v[svert[2]] + s_2*v[svert[3]];
		}
	}

	return average;

}

Mesh::Vertex getDivisionPoint(const Mesh &m, size_t inda, size_t indb)
{
	Mesh::Vertex av1 = getAveragedVertex(m, inda, indb);
	Mesh::Vertex av2 = getAveragedVertex(m, indb, inda);

	Mesh::Vertices v = m.get_all_vertices();

	if(av1 == Mesh::Vertex(0,0,0) || av2 == Mesh::Vertex(0,0,0))
		return (v[inda]+v[indb])/2;

	return getAveragedVertex(m, inda, indb) + getAveragedVertex(m, indb, inda);
}


Mesh surfaces::mButterflySubdivision(const Mesh &source, int iterations)
{	
	Mesh msrc = source;

	for(int t = 0; t < iterations; ++t)
	{
		Mesh mdst = Mesh(msrc.get_all_vertices().size());
		
		//Maps edges (indices of two points) with their division points
		std::map<std::pair<size_t, size_t>, size_t> division_map;
		
		Mesh::Vertices vertices = msrc.get_all_vertices();

		//Add the old vertices into the new mesh
		for(size_t ind = 0; ind < vertices.size(); ++ind)
		{			
			mdst.add_vertex(vertices[ind]);
		}

		//Add new subdivision vertices into the new mesh
		for(size_t ind = 0; ind < vertices.size(); ++ind)
		{			
			Mesh::AdjacentVerticesPerVertex adjacent = msrc.get_neighbouring_vertices(ind);	
			Mesh::AdjacentVerticesPerVertex::const_iterator it_adjacent = adjacent.begin();
			while(it_adjacent != adjacent.end())
			{
				size_t ind_adjacent = *it_adjacent;
												
				//"Id" for the edge [ind, ind_adjacent], smaller index first 
				std::pair<size_t, size_t> keyEdge(ind < ind_adjacent ? ind : ind_adjacent, 
												  ind > ind_adjacent ? ind : ind_adjacent);
								
				//Check: if edge [ind, ind_adjacent] is not yet processed.
				if(division_map.find(keyEdge) == division_map.end())
				{
					//Calculate a division vertex for [ind, ind_adjacent].
					Mesh::Vertex vertex_div = getDivisionPoint(msrc, ind, ind_adjacent);
					//Add the division vertex into the new mesh
					size_t ind_div = mdst.add_vertex(vertex_div);
					//Map the edge [ind, ind_adjacent] with the division vertex
					division_map[keyEdge]=ind_div;
				}

				++it_adjacent;
			}
						
		}

		Mesh::Faces faces = msrc.get_all_faces();

		//Create faces of the new mesh
		Mesh::Faces::const_iterator it = faces.begin();
		while(it != faces.end())
		{ 
			
			//Get indices of the division vertices
			size_t inds_div[3];
			for(int c = 0; c < 3; ++c)
			{
				size_t p1 = 0, p2 = 0;
				
				//Get the current edge of the face
				switch(c)
				{
					case 0: p1 = it->A();
							p2 = it->B();
							break;
					case 1: p1 = it->B();
							p2 = it->C();
							break;
					case 2: p1 = it->C();
							p2 = it->A();
				}

				//Generate the "Id" of the current edge
				std::pair<size_t, size_t> keyEdge(p1 < p2 ? p1 : p2,
                                                  p1 > p2 ? p1 : p2);
				//Get the index of the corresponding division vertex
				inds_div[c] = division_map[keyEdge];
			}
			
			//Generate new four faces.
			mdst.add_face(Mesh::Face(it->A(), inds_div[2], inds_div[0]));
			mdst.add_face(Mesh::Face(it->B(), inds_div[0], inds_div[1]));
			mdst.add_face(Mesh::Face(it->C(), inds_div[1], inds_div[2]));
			mdst.add_face(Mesh::Face(inds_div[0], inds_div[1], inds_div[2]));

			++it;
		}

		msrc = mdst;
	}

	return msrc;
}

} //namespace methods
} //namespace common
