// Using this file since the original fempic mesh loaders and functions are used

#include "fempic.h"
#include <chrono>
#include <ctime>

int Particles::get_num_particles()
{ 
	return particle_cell_index.size(); 
}

void Particles::collect_particle_to_add(
	const int cell_index, 
	array<double,DIMENSIONS> pos, 
	array<double,DIMENSIONS> vel, 
	array<double,DIMENSIONS> ef, 
	array<double,NODES_PER_CELL> lc)
{
	particle_cell_index_to_add.push_back(cell_index);
	particle_pos_to_add.push_back(pos);	
	particle_vel_to_add.push_back(vel);
	particle_ef_to_add.push_back(ef);
	particle_lc_to_add.push_back(lc);
}

void Particles::insert_collected_particles()
{
	for (int i = 0; i < (int)particle_cell_index_to_add.size(); i++)
	{
		particle_cell_index.push_back(particle_cell_index_to_add[i]);
		particle_pos.push_back(particle_pos_to_add[i]);	
		particle_vel.push_back(particle_vel_to_add[i]);
		particle_ef.push_back(particle_ef_to_add[i]);
		particle_lc.push_back(particle_lc_to_add[i]);
	}
	printf("insert_collected_particles Allocating set [PARTICLES] with size [%d]\n", (int)particle_cell_index_to_add.size());
	particle_cell_index_to_add.clear();
	particle_pos_to_add.clear();
	particle_vel_to_add.clear();
	particle_ef_to_add.clear();
	particle_lc_to_add.clear();
}

void Particles::mark_to_remove_particle(const int particle_index)
{
	particle_indexes_to_remove.push_back(particle_index);
}

void Particles::remove_marked_particles()
{
	int remove_count = 0;
	for (auto& particle_index : particle_indexes_to_remove)
	{
		particle_pos.erase(particle_pos.begin() + particle_index - remove_count);
		particle_vel.erase(particle_vel.begin() + particle_index - remove_count);
		particle_ef.erase(particle_ef.begin() + particle_index - remove_count);
		particle_lc.erase(particle_lc.begin() + particle_index - remove_count);
		particle_cell_index.erase(particle_cell_index.begin() + particle_index - remove_count);
		
		remove_count++;
	}
	printf("remove_marked_particles set [PARTICLES] with size [%d]\n", (int)particle_indexes_to_remove.size());
	particle_indexes_to_remove.clear();
}


/*loads and initializes volume mesh*/
bool LoadVolumeMesh(const string file_name, Volume &volume)
{
	/*open file*/
	ifstream in(file_name);
	if (!in.is_open()) {cerr<<"LOAD : Failed to open "<<file_name<<endl; return false;}
	
	/*read number of nodes and cells*/
	int n_nodes, n_cells;
	in>>n_nodes>>n_cells;
	cout<<"LOAD : Mesh contains "<<n_nodes<<" nodes and "<<n_cells<<" cells"<<endl;

	/*read the nodes*/
	for (int n=0;n<n_nodes;n++)
	{
		int index;
		double x, y, z;

		in >> index >> x >> y >> z;
		if (index!=n+1) cout<<"LOAD : Inconsistent node numbering"<<endl;
			
		volume.nodes.emplace_back(x/1000.,y/1000.,z/1000.);		
	}
	
	/*read cells, this will also contain edges and triangles*/
	for (int e=0;e<n_cells;e++)
	{
		int index, type;
		int n1, n2, n3, n4;

		in >> index >> type;
		
		if (type!=304) {string s; getline(in,s);continue;} // after line 4775 in mesh.dat
		
		in >> n1 >> n2 >> n3 >> n4;

		/*flipping nodes 2 & 3 to get positive volumes*/
		volume.cells.emplace_back(n1-1, n2-1, n3-1, n4-1);		
	}

	/*reset number of nodes and cells since we skipped bunch of lines and triangles*/
	n_nodes = volume.nodes.size();
	n_cells = volume.cells.size();

	cout<<"LOAD : Final Mesh volume "<<n_nodes<<" nodes and "<<n_cells<<" cells"<<endl;

	/*compute cell volumes*/
	for (Tetra &tet: volume.cells)
	{
		double M[4][4];
		
		/*set first column to 1*/
		for (int i=0;i<4;i++) M[i][0] = 1;
		
		/*loop over vertices, i.e nodes???*/
		for (int v=0;v<4;v++)
		{
			for (int dim=0;dim<3;dim++)
			{
				M[0][dim+1] = volume.nodes[tet.con[0]].pos[dim];
				M[1][dim+1] = volume.nodes[tet.con[1]].pos[dim];
				M[2][dim+1] = volume.nodes[tet.con[2]].pos[dim];
				M[3][dim+1] = volume.nodes[tet.con[3]].pos[dim];		
			}
		}
		
		/*volume is (1/6)*det4(M)*/
		tet.volume = (1.0/6.0)*det4(M);
		
		/*flip ABCD to ADBC if negative volume*/
		if (tet.volume < 0) 
		{
			int t=tet.con[1];
			tet.con[1]=tet.con[3];
			tet.con[3]=t;
			tet.volume=-tet.volume;
		}
	}
	
	/*precompute 3x3 determinants for LC computation*/
	for (Tetra &tet:volume.cells)
	{
		double M[3][3];
		/*loop over vertices*/
		for (int v=0; v<4; v++)
		{
			int v2,v3,v4;
			
			switch (v)
			{
				case 0: v2=1; v3=2; v4=3; break;
				case 1: v2=3; v3=2; v4=0; break;
				case 2: v2=3; v3=0; v4=1; break;
				case 3: v2=1; v3=0; v4=2; break;				
			}
			
			double *p2 = volume.nodes[tet.con[v2]].pos;
			double *p3 = volume.nodes[tet.con[v3]].pos;
			double *p4 = volume.nodes[tet.con[v4]].pos;
			
			/*alpha*/
			M[0][0] = p2[0];
			M[0][1] = p2[1];
			M[0][2] = p2[2];
			M[1][0] = p3[0];
			M[1][1] = p3[1];
			M[1][2] = p3[2];
			M[2][0] = p4[0];
			M[2][1] = p4[1];
			M[2][2] = p4[2];
			tet.alpha[v] = det3(M);
			
			/*beta*/
			M[0][0] =1;
			M[0][1] = p2[1];
			M[0][2] = p2[2];
			M[1][0] = 1;
			M[1][1] = p3[1];
			M[1][2] = p3[2];
			M[2][0] = 1;
			M[2][1] = p4[1];
			M[2][2] = p4[2];
			tet.beta[v] = det3(M);

			/*gamma*/
			M[0][0] =1;
			M[0][1] = p2[0];
			M[0][2] = p2[2];
			M[1][0] = 1;
			M[1][1] = p3[0];
			M[1][2] = p3[2];
			M[2][0] = 1;
			M[2][1] = p4[0];
			M[2][2] = p4[2];
			tet.gamma[v] = det3(M);
			
			/*delta*/
			M[0][0] =1;
			M[0][1] = p2[0];
			M[0][2] = p2[1];
			M[1][0] = 1;
			M[1][1] = p3[0];
			M[1][2] = p3[1];
			M[2][0] = 1;
			M[2][1] = p4[0];
			M[2][2] = p4[1];
			tet.delta[v] = det3(M);
		}
	}

	/*build cell connectivity, there is probably a faster way*/
	cout<<"LOAD : Building cell connectivity"<<endl;
	
	/*reset connectivities*/
	for (int l = 0; l < n_cells; l++)
	{
		Tetra &tet = volume.cells[l];
		for (int v=0; v<4; v++) 
		{
			tet.cell_con[v] = -1;	/*no neighbor*/		
		}	
	}


	for (int l = 0; l < n_cells; l++)
	{
		Tetra &tet = volume.cells[l];
		int v1,v2,v3;
		for (int v=0; v<4; v++)
		{
			/*skip if already set*/
			if (tet.cell_con[v] >= 0) continue;

			switch(v)
			{
				case 0: v1=1; v2=2; v3=3; break;
				case 1: v1=2; v2=3; v3=0; break;
				case 2: v1=3; v2=0; v3=1; break;
				case 3: v1=0; v2=1; v3=2; break;
			}
			
			/*loop over the tets again looking for one with these three vertices*/
			for (int m=l+1; m<n_cells; m++)
			{
				Tetra &other = volume.cells[m];
				
				bool matches[4] = {false,false,false,false};
				int count = 0;
				for (int k=0;k<4;k++)
				{
					if (other.con[k]==tet.con[v1] ||
						other.con[k]==tet.con[v2] ||
						other.con[k]==tet.con[v3]) 
						{
							count++;matches[k]=true;
						}
				}
				
				/*if three vertices match*/
				if (count==3) 
				{
					tet.cell_con[v] = m;		

					/*set the cell connectivity for the index without a matching vertex to l, l is itself*/
					for (int k=0;k<4;k++)
						if(!matches[k]) other.cell_con[k] = l;
				}
			}
		}		
	}

	/*also compute node volumes by scattering cell volumes,this can only be done after 3x3 dets are computed*/
	
	/*first set all to zero*/
	for (Node &node : volume.nodes) 
	{
		node.volume = 0;
	}

	for (int i=0; i<n_cells; i++)
	{
		Tetra &tet = volume.cells[i];
		
		int part_cell_index = i;
		std::array<double,DIMENSIONS> part_pos 		= {0.0,0.0,0.0};
		std::array<double,NODES_PER_CELL> part_lc 	= {0.0,0.0,0.0,0.0};

		/*compute centroid position*/
		for (int dim=0;dim<DIMENSIONS;dim++)
		{
			for (int v=0;v<4;v++)
			{
				part_pos[dim] += 0.25*volume.nodes[tet.con[v]].pos[dim];	
			}		
		}
		
		bool found = XtoLtet(part_cell_index,part_lc,part_pos,volume,false);
		if (!found) cout<<"LOAD : something is wrong"<<endl;
		
		for (int v=0; v<4; v++)
		{
			volume.nodes[tet.con[v]].volume += part_lc[v]*tet.volume;
		}
		
	}

	/*mark nodes on open faces as open*/
	for (size_t e=0; e<volume.cells.size(); e++)
	{
		Tetra &tet = volume.cells[e];
		for (int v=0;v<4;v++)
		{
			if (tet.cell_con[v]<0)	/*no neighbor*/
			{
				for (int i=0;i<4;i++)
				{
					if (i!=v) volume.nodes[tet.con[i]].type=OPEN;
				}
			}
		}
	}
	
	cout<<"LOAD : Done loading "<<file_name<<endl;	
	return true;
}

/*loads nodes from a surface mesh file and sets them to the specified node type*/
bool LoadSurfaceMesh(const string file_name, Volume &volume, NodeType node_type)
{
	/*open file*/
	ifstream in(file_name);
	if (!in.is_open()) {cerr<<"LOAD : Failed to open "<<file_name<<endl; return false;}
	
	/*read number of nodes and cells*/
	int n_nodes, n_cells;
	in>>n_nodes>>n_cells;
	cout<<"LOAD : Mesh contains "<<n_nodes<<" nodes and "<<n_cells<<" cells"<<endl;

	int nn = volume.nodes.size();

	/*read the nodes*/
	for (int n=0;n<n_nodes;n++)
	{
		int index;
		double x, y, z;

		in >> index >> x >> y >> z;
		
		if (index<1 || index>nn) {cerr<<"LOAD : Incorrect node number "<<index<<endl;continue;}
		volume.nodes[index-1].type=node_type;
	}
		
	cout<<"LOAD : Done loading "<<file_name<<endl;	
	return true;
}

//*************************************************************************************************
/*converts physical coordinate to logical
Returns true if particle matched to a tet
*/
bool XtoLtet(
	int& part_cell_index, 
	array<double, NEIGHBOUR_CELLS>& part_lc, 
	const array<double, DIMENSIONS>& part_pos, 
	const Volume &volume, 
	bool search)
{	
	/*first try the current tetrahedron*/
	const Tetra &tet = volume.cells[part_cell_index];
	
	bool inside = true;
	
	for (int i=0; i<NEIGHBOUR_CELLS; i++) /*loop over vertices*/
	{
		part_lc[i] = (1.0/6.0)*(tet.alpha[i] - part_pos[0]*tet.beta[i] + 
					  part_pos[1]*tet.gamma[i] - part_pos[2]*tet.delta[i])/tet.volume;
		
		if (part_lc[i]<0 || part_lc[i]>1.0) 
			inside=false;
	}	
	
	if (inside) return true;
	
	if (!search) return false;

	/*we are outside the last known tet, find most negative weight*/
	int min_i=0;
	double min_lc=part_lc[0];
	
	for (int i=1; i<NEIGHBOUR_CELLS; i++)
	{
		if (part_lc[i]<min_lc) 
		{
			min_lc=part_lc[i];
			min_i=i;
		}
	}
	
	/*is there a neighbor in this direction?*/
	if (tet.cell_con[min_i]>=0)
	{
		part_cell_index = tet.cell_con[min_i];
		return XtoLtet(part_cell_index, part_lc, part_pos, volume);
	}
	
	return false;
}

/*** FUNCTIONS ***/

/*saves volume mesh*/
void OutputMesh(int ts, Volume &volume, double *phi, double *ef, double *ion_den)
{
	stringstream ss;
	ss<<"results/mesh_"<<setfill('0')<<setw(4)<<ts+1<<".vtu";
	ofstream out(ss.str());
	if (!out.is_open()) {cerr<<"Failed to open file "<<ss.str()<<endl;exit(-1);}
	
	/*header*/
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out<<"<UnstructuredGrid>\n";
	out<<"<Piece NumberOfPoints=\""<<volume.nodes.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
	out<<"NumberOfStrips=\"0\" NumberOfCells=\""<<volume.cells.size()<<"\">\n";
	
	/*points*/
	out<<"<Points>\n";
	out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (Node &node: volume.nodes)
		out<<std::fixed << setprecision(PRINT_PRECISION)<<node.pos[0]<<" "<<node.pos[1]<<" "<<node.pos[2]<<"\n";		
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	/*Cells*/
	out<<"<Cells>\n";
	out<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (Tetra &tetra: volume.cells)
		out<<tetra.con[0]<<" "<<tetra.con[1]<<" "<<tetra.con[2]<<" "<<tetra.con[3]<<"\n";
	out<<"</DataArray>\n";
	
	out<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (size_t e=0; e<volume.cells.size();e++)
		out<<(e+1)*4<<" ";
	out<<"\n";
	out<<"</DataArray>\n";
	
	out<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (size_t e=0; e<volume.cells.size();e++)
		out<<"10 ";
	out<<"\n";
	out<<"</DataArray>\n";		
	out<<"</Cells>\n";

	/*save point data*/
	out<<"<PointData Scalars=\"phi\">\n";
	out<<"<DataArray type=\"Int32\" Name=\"node_index\" format=\"ascii\">\n";
	for (size_t n=0; n<volume.nodes.size();n++)
		out<<n<<" ";
	out<<"\n";
	out<<"</DataArray>\n";
	out<<"<DataArray type=\"Int32\" Name=\"node_type\" format=\"ascii\">\n";
	for (size_t n=0; n<volume.nodes.size();n++)
		out<<volume.nodes[n].type<<" ";
	out<<"\n";
	out<<"</DataArray>\n";
	
	out<<"<DataArray type=\"Float32\" Name=\"phi\" format=\"ascii\">\n";
	for (size_t n=0; n<volume.nodes.size();n++)
		out<<std::fixed << setprecision(PRINT_PRECISION)<<phi[n]<<" ";
	out<<"\n";
	out<<"</DataArray>\n";
	
	out<<"<DataArray type=\"Float32\" Name=\"ion_den\" format=\"ascii\">\n";
	for (size_t n=0; n<volume.nodes.size();n++)
		out<<std::fixed << setprecision(PRINT_PRECISION)<<ion_den[n]<<" ";
	out<<"\n";
	out<<"</DataArray>\n";

	out<<"</PointData>\n";

	/*save cell data*/
	out<<"<CellData Vectors=\"ef\">\n";
	out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"ef\" format=\"ascii\">\n";
	for (size_t e=0; e<volume.cells.size();e++)
		out<<std::fixed << setprecision(PRINT_PRECISION)<<ef[e * DIMENSIONS + 0]<<" "<<ef[e * DIMENSIONS + 1]<<" "<<ef[e * DIMENSIONS + 2]<<" ";
	out<<"\n";
	out<<"</DataArray>\n";
	
	out<<"<DataArray type=\"Float32\" Name=\"cell_volume\" format=\"ascii\">\n";
	for (Tetra &tet:volume.cells)
		out<<tet.volume<<" ";
	out<<"\n";
	out<<"</DataArray>\n";

	out<<"</CellData>\n";

	out<<"</Piece>\n";
	out<<"</UnstructuredGrid>\n";
	out<<"</VTKFile>\n";

	out.close();
}


/*saves particle data*/
void OutputParticles(int ts, Particles &particles)
{
	stringstream ss;
	ss<<"results/particles_"<<setfill('0')<<setw(4)<<ts+1<<".vtp";
	ofstream out(ss.str());

	if (!out.is_open()) {cerr<<"Failed to open output file "<<endl;exit(-1);}
	
	/*header*/
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out<<"<PolyData>\n";
	out<<"<Piece NumberOfPoints=\""<<particles.particle_cell_index.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
	out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";
	
	/*points*/
	out<<"<Points>\n";
	out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (auto &part: particles.particle_pos)
		out<<std::fixed << setprecision(PRINT_PRECISION)<<part[0]<<" "<<part[1]<<" "<<part[2]<<"\n";		
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	out<<"</Piece>\n";
	out<<"</PolyData>\n";
	out<<"</VTKFile>\n";

	out.close();	
}

/*computes determinant of a 4x4 matrix*/
double det4(double (*M)[4])
{
	double M0[3][3];
	double M1[3][3];
	double M2[3][3];
	double M3[3][3];
	
	for (int i=0;i<3;i++)
	{
		M0[i][0]=M[i+1][1];
		M0[i][1]=M[i+1][2];
		M0[i][2]=M[i+1][3];
		
		M1[i][0]=M[i+1][0];
		M1[i][1]=M[i+1][2];
		M1[i][2]=M[i+1][3];

		M2[i][0]=M[i+1][0];
		M2[i][1]=M[i+1][1];
		M2[i][2]=M[i+1][3];
		
		M3[i][0]=M[i+1][0];
		M3[i][1]=M[i+1][1];
		M3[i][2]=M[i+1][2];
	}
	
	return M[0][0]*det3(M0) - 
		   M[0][1]*det3(M1) +
		   M[0][2]*det3(M2) -
		   M[0][3]*det3(M3);
}

/*computes determinant of a 3x3 matrix*/
double det3(double (*M)[3])
{
	return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])- 
		   M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+
		   M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);	
}

/*helper functions for matrix math, y=A*x */
void matVecMultiply(double *y, double**A, double *x, int nu)
{
	for (int i=0; i<nu; i++)
	{
		y[i] = 0;
		for (int j=0; j<nu; j++)
			y[i] += A[i][j]*x[j];
	}
}

/*computes y=v1-v2*/
void vecVecSubtract(double *y, double *v1, double *v2, int nu)
{
	for (int i=0;i<nu;i++)
			y[i] = v1[i]-v2[i];
}

//*************************************************************************************************
void EnrichArrays(
	Volume& volume, 
	double *&node_bnd_pot_tmp, double *&node_pot_tmp, double *&node_pos_tmp, double *&node_volume_tmp, double *&node_ion_den_tmp,
	int *&cell_to_nodes_tmp, int *&cell_to_cell_tmp, double *&cell_det_tmp, double *&cell_volume_tmp, double *&cell_ef_tmp)
{

	if (!LoadVolumeMesh("dat_files/mesh.data",volume) ||
		!LoadSurfaceMesh("dat_files/inlet.data",volume,INLET) ||
		!LoadSurfaceMesh("dat_files/sphere.data",volume,SPHERE))
		{
			return;
		}

	int n_nodes 	= volume.nodes.size();
	int n_cells 	= volume.cells.size();

	cell_ef_tmp 		= new double[n_cells * DIMENSIONS];
	node_bnd_pot_tmp 	= new double[n_nodes];
	node_pot_tmp 		= new double[n_nodes];
	node_ion_den_tmp 	= new double[n_nodes];

	cell_to_nodes_tmp 	= new int[n_cells * NODES_PER_CELL];
	cell_to_cell_tmp	= new int[n_cells * NEIGHBOUR_CELLS];
	cell_det_tmp		= new double[n_cells * DET_FIELDS * NEIGHBOUR_CELLS]; // arranged as [alpha,beta,gamma,delta] * 4 neighbours
	cell_volume_tmp		= new double[n_cells];

	node_pos_tmp 		= new double[n_nodes * DIMENSIONS];
	node_volume_tmp		= new double[n_nodes];

	for (int n=0; n<n_nodes; n++)
	{
		switch (volume.nodes[n].type)
		{
			case INLET: 	node_bnd_pot_tmp[n] = 0; break;		/*phi_inlet*/
			case SPHERE: 	node_bnd_pot_tmp[n] = -100; break; 	/*phi_sphere*/
			default: 		node_bnd_pot_tmp[n] = 0;			/*default*/
		}

		node_pot_tmp[n] = 0.0f;

		Node &node = volume.nodes[n];
		
		for (int dim=0; dim<DIMENSIONS; dim++)
			node_pos_tmp[n * DIMENSIONS + dim] = node.pos[dim];
	
		node_volume_tmp[n] 	= node.volume;
	}

	for (int cellID=0; cellID<n_cells; cellID++)
	{
		Tetra &tet = volume.cells[cellID];
		
		for (int nodeCon=0; nodeCon<NODES_PER_CELL; nodeCon++)
			cell_to_nodes_tmp[cellID * NODES_PER_CELL + nodeCon] = tet.con[nodeCon];
		
		for (int cellCon=0; cellCon<NEIGHBOUR_CELLS; cellCon++)
		{
			cell_to_cell_tmp[cellID * NEIGHBOUR_CELLS + cellCon] 	= tet.cell_con[cellCon];

			cell_det_tmp[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 0] = tet.alpha[cellCon];
			cell_det_tmp[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 1] = tet.beta[cellCon];
			cell_det_tmp[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 2] = tet.gamma[cellCon];
			cell_det_tmp[(cellID * NEIGHBOUR_CELLS + cellCon) * DET_FIELDS + 3] = tet.delta[cellCon];
		}

		cell_volume_tmp[cellID] = tet.volume;

		cell_ef_tmp[cellID * DIMENSIONS + 0] = 0.0;
		cell_ef_tmp[cellID * DIMENSIONS + 1] = 0.0;
		cell_ef_tmp[cellID * DIMENSIONS + 2] = 0.0;
	}
}

//*************************************************************************************************
void print_particles_m(const double* part_pos, const double* part_vel, const int* cell_index, int n_particles, const std::string suffix)
{
	printf("PRINT : print_particles_m %s ... n_particles %d\n", suffix.c_str(), n_particles);

	std::string file_name1 = std::string("files/") + suffix + "_particle_position.dat"; 
	std::string file_name2 = std::string("files/") + suffix + "_particle_velocity.dat"; 
	std::string file_name3 = std::string("files/") + suffix + "_particle_other.dat"; 

	FILE *fp[3];
	fp[0] = fopen(file_name1.c_str(), "w");
	fp[1] = fopen(file_name2.c_str(), "w");
	fp[2] = fopen(file_name3.c_str(), "w");

	for (int i = 0; i < n_particles; i++)
	{
		fprintf(fp[0], "%d,%+2.30lE,%+2.30lE,%+2.30lE\n", i, part_pos[i * DIMENSIONS + 0], part_pos[i * DIMENSIONS + 1], part_pos[i * DIMENSIONS + 2]);
		fprintf(fp[1], "%d,%+2.30lE,%+2.30lE,%+2.30lE\n", i, part_vel[i * DIMENSIONS + 0], part_vel[i * DIMENSIONS + 1], part_vel[i * DIMENSIONS + 2]);
		fprintf(fp[2], "%d,%d\n", i, cell_index[i]);
	}

	fclose(fp[0]);
	fclose(fp[1]);
	fclose(fp[2]);
}

//*************************************************************************************************
void print_fields_m(double *ion_den, double *field_potential, double *electric_field, int n_nodes, int n_cells, const std::string suffix)
{
	printf("PRINT : print_fields_m %s ... num_nodes %d | num cells %d\n", suffix.c_str(), n_nodes, n_cells);

	std::string file_name1 = std::string("files/") + suffix + "_field_pot_den.dat"; 
	std::string file_name2 = std::string("files/") + suffix + "_field_ef.dat"; 

	FILE *fp[2];
	fp[0] = fopen(file_name1.c_str(), "w");
	fp[1] = fopen(file_name2.c_str(), "w");

	for (int i = 0; i < n_nodes; i++)
	{
		fprintf(fp[0], "%d,%+2.30lE,%+2.30lE\n", i, ion_den[i], field_potential[i]);
	}

	for (int i = 0; i < n_cells; i++)
	{
		// 	printf("printing file - %d\t%+2.30lE\t%+2.30lE\t%+2.30lE\n\n", i, electric_field[i * DIMENSIONS + 0], electric_field[i * DIMENSIONS + 1],electric_field[i * DIMENSIONS + 2]);		
		fprintf(fp[1], "%d,%+2.30lE,%+2.30lE,%+2.30lE\n", i, electric_field[i*DIMENSIONS], electric_field[i*DIMENSIONS+1], electric_field[i*DIMENSIONS+2]);
	}
	
	fclose(fp[0]);
	fclose(fp[1]);
}

//*************************************************************************************************
std::mt19937 mt_gen(0);		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd() {
	return rnd_dist(mt_gen);
}

//*************************************************************************************************
int get_num_particles_to_inject(double dt, double& remainder)
{
	double fnum_mp 				= ((num_per_sec * dt) / OP_CONST_spwt + remainder);	/*fraction number of macroparticles*/  
	int num_particles_to_add 	= (int)fnum_mp;										/*integer number of macroparticles*/
	remainder 					= (fnum_mp - num_particles_to_add); 				/*update reminder*/

	return num_particles_to_add;
}

//*************************************************************************************************

