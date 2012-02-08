#include "common/methods/mbutterfly_mesh_subdiv.hpp"
#include <map>
#include <list>

namespace common
{

namespace methods
{

//Returns a consequence of vertices that defines a stencil (clockwise or 
//counter-clockwise traversal)
std::vector<size_t> getRelativeStencil(const Mesh &m, size_t inda, size_t indb)
{
    //Array of the sorted adjacent vertices that define a stencil
    std::vector<size_t> svert;

	Mesh::Faces f = m.get_all_faces();
	
	Mesh::AdjacentFacesPerVertex adjacentFaces = m.get_neighbouring_faces_by_vertex(inda);
	Mesh::AdjacentVerticesPerVertex adjacentVertices = m.get_neighbouring_vertices(inda);

	//If the central vertex is on the boundary, calculate
	//the average value using the boundary schema
	if(adjacentFaces.size() != adjacentVertices.size())
	{
		return svert;
	}
		
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

    return svert;
}

//Returns the average value of vertices for irregular K-stencils
Mesh::Vertex getKStencilAverage(const Mesh &m, std::vector<size_t> sten, Mesh::Vertex center)
{
    Mesh::Vertices v = m.get_all_vertices();
    Mesh::Vertex average(0,0,0);

    size_t k = sten.size();
   
    //These cases are described in the reference paper (see .hpp).
    if (k >= 5)
    {
        float pi = 3.141592f;
        for(size_t j = 0; j < k; ++j)
        {
            float s_j = (0.25f + cos(2*pi*j/k) + 0.5f*cos(4*pi*j/k))/k;
            average += s_j*(v[sten[j]]-center);
        }
    }
    else if (k == 4)
    {
        float s_0, s_1, s_2, s_3;
        s_0 = 3/8.0f;
        s_1 = 0;
        s_2 = -1/8.0f;
        s_3 = 0;

        average = s_0*(v[sten[0]]-center) + s_1*(v[sten[1]]-center) +
                  s_2*(v[sten[2]]-center) + s_3*(v[sten[3]]-center);
    }
    else if(k == 3)
    {
        float s_0, s_1, s_2;
        s_0 = 5/12.0f;
        s_1 = -1/12.0f;
        s_2 = s_1;

        average = s_0*(v[sten[0]]-center) + s_1*(v[sten[1]]-center) +
                  s_2*(v[sten[2]]-center);
    }

    average += center;

    return average;
}

//Returns the point that divides the given edge [inda, indb] 
Mesh::Vertex getDivisionPoint(const Mesh &m, size_t inda, size_t indb)
{
    Mesh::Vertices v = m.get_all_vertices();

    //Calculate two stencils for the end points of the edge
	std::vector<size_t> sten1 = getRelativeStencil(m, inda, indb);
	std::vector<size_t> sten2 = getRelativeStencil(m, indb, inda);

    //Get the valence of the stencils' centers
	int k1 = sten1.size();
    int k2 = sten2.size();
    
    Mesh::Vertex average(0,0,0);

	//The Modified Butterfly subdivision surface scheme implementation
    //(see the reference paper, .hpp)
    {        
        bool isOrdered = false;

        //Border case: two border vertices. TODO: need to be improved
        if(k1 == 0 && k2 == 0)
        {
            average = (v[inda]+v[indb])/2;
        }
        //Border case: border and inner vertices
        else if((isOrdered = (k1 == 0 && k2 != 0)) == true || (k1 != 0 && k2 == 0)) //a trick with C4706 warning
        {
            size_t indc = indb;

            if(!isOrdered)
            {
                std::vector<size_t> tmp = sten2;
                sten2 = sten1;
                sten1 = tmp;
                indc = inda;
            }

            average = getKStencilAverage(m, sten2, v[indc]);
        }
	    //Regular case. Two vertices of valence 6
        else if(k1 == 6 && k2 == 6) 
	    {
            //Regular stencil weights
            float a, b, c, d;
            a = 1/2.0f;
            b = 1/8.0f;
            c = -1/16.0f;
            d = 0.0f;

		    average = a*(v[inda] + v[indb]) +
                      b*(v[sten1[1]] + v[sten1[5]]) +
                      c*(v[sten1[2]] + v[sten1[4]] + v[sten2[2]] + v[sten2[4]]) +
                      d*(v[sten1[3]] + v[sten2[3]]);
	    }
        //Non-regular cases (modification of the Butterfly scheme):
        //a) 6-vertex and K-vertex, K!=6 
	    else if((isOrdered = (k1 == 6 && k2 != 6)) == true || (k1 != 6 && k2 == 6)) //a trick with C4706 warning
	    {
            if(!isOrdered)
            {
                std::vector<size_t> tmp = sten2;
                sten2 = sten1;
                sten1 = tmp;
            }

            average = getKStencilAverage(m, sten2, v[sten1[0]]);
	    }
        //b) K-vertex and K-vertex, K!=6
        else
        {
            average = (getKStencilAverage(m, sten2, v[sten1[0]]) + 
                       getKStencilAverage(m, sten1, v[sten2[0]])) / 2;
        }
    }

    return average;

}

//Returns a subdivided mesh acquired from the source mesh by application of
//the given number of subdivision iterations using the Modified Butterfly 
//scheme
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
