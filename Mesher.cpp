//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Mesher.cpp
//OBJECTIVE:	Generate the mesh of RVE, including cracks.
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include "Mesher.h"

//---------------------------------------------------------------------------
//�������ɺ���
int Mesher::Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, const string disc, const struct Weight_func &weight_func)
{
	//���ɱ�������������
	if(Brick_background_mesh(grid_size, geom_rve, disc, weight_func)== 0) return 0;
	//��������Ԫ��������
	if(Element_material_property(ele_prop)==0) return 0;
	//����������������
	if(Add_crack_to_element(Cracks, geom_rve)==0) return 0;
	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Generate meshes for damage and fracture propagation
int Mesher::Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Geom_RVE &geom_rve, const string disc, const vector<bool> &full_dam_eles)
{
	//Generate brick grids
	if(Brick_background_mesh(grid_size, geom_rve, disc, full_dam_eles)== 0) return 0;
	//��������Ԫ��������
	if(Element_material_property(ele_prop)==0) return 0;
	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Generate meshes with cracks for damage and fracture propagation
int Mesher::Generate_mesh(const struct Grid_size &grid_size, const struct Ele_prop &ele_prop, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, const string disc, const vector<bool> &full_dam_eles)
{
	//���ɱ�������������
	if(Brick_background_mesh(grid_size, geom_rve, disc, full_dam_eles)== 0) return 0;
	//��������Ԫ��������
	if(Element_material_property(ele_prop)==0) return 0;
	//����������������
	if(Add_crack_to_element(Cracks, geom_rve)==0) return 0;
	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh data(DGFEM)
int Mesher::Import_mesh_reconfiguration(string output_file_name, const string disc, const vector<bool> &full_dam_eles)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Read_mesh(output_file_name)== 0) return 0;

	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Read_mesh(output_file_name)== 0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		if(full_dam_eles.size()!=0)
		{
			//------------------------------------------
			Deter_nodes_relative_eles();

			//------------------------------------------
			vector<Node> nodesFEM(nodes);
			nodes.clear();
			vector<int> new_count;
			//------------------------------------------
			//Define the DG nodes
			for(int i=0; i<(int)nodesFEM.size(); i++)
			{
				//-----------------------------------------------------------------------
				//Damaged node, at least one of its relative elements is damaged.
				//for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
				//{
				//	if(full_dam_eles[nodesFEM[i].relative_eles[j]]) 
				//	{
				//		nodesFEM[i].DG_node = true;
				//		break;
				//	}
				//}

				//-----------------------------------------------------------------------
				//Full damaged node, Its relative elements are all damaged. 
				nodesFEM[i].DG_node = true;
				for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
				{
					if(!full_dam_eles[nodesFEM[i].relative_eles[j]]) 
					{
						nodesFEM[i].DG_node = false;
						break;
					}
				}

				if(nodesFEM[i].DG_node) new_count.push_back(-1);
				else
				{
					new_count.push_back((int)nodes.size());
					nodes.push_back(nodesFEM[i]);
				}
			}

			//------------------------------------------
			for(int i=0; i<(int)elements.size(); i++)
			{
				elements[i].flag = 1;
				elements[i].DG_elem=false;
				int counter=0;
				for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
				{
					if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
					{
						counter++;
						nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
						elements[i].nodes_id[j] = (int)nodes.size()-1;
					}
					else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
					{
						int pp = elements[i].nodes_id[j];
						elements[i].nodes_id[j] = new_count[pp];
					}
					else
					{
						cout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;	
						hout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;
						return 0;
					}
				}
				if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
			}

			//--------------------------------------------
			//clear the relative elements of nodes
			for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_eles.clear();
		}
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh data(DGFEM) for pure nonlocal or pure nonlocal fracture propagation
int Mesher::Import_mesh_reconfiguration(string output_file_name, const string disc, const struct Weight_func &weight_func)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Read_mesh(output_file_name)== 0) return 0;

	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Read_mesh(output_file_name)== 0) return 0;

		//------------------------------------------
		vector<Node> nodesFEM(nodes);
		nodes.clear();
		vector<int> new_count;
		WeightFunc Wei_Fun;  //the class of weighting function
		//------------------------------------------
		//Define the DG nodes
		for(int i=0; i<(int)nodesFEM.size(); i++)
		{
			Point_3D poi_tem(nodesFEM[i].x, nodesFEM[i].y, nodesFEM[i].z);
			if(fabs(Wei_Fun.Value_weight_function(poi_tem, weight_func)-1.0)<=Zero)
			{
				nodesFEM[i].DG_node = true;
				new_count.push_back(-1);
			}
			else
			{
				new_count.push_back((int)nodes.size());
				nodes.push_back(nodesFEM[i]);
			}
		}

		//------------------------------------------
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].flag = 1;
			elements[i].DG_elem=false;
			int counter=0;
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
				{
					counter++;
					nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
					elements[i].nodes_id[j] = (int)nodes.size()-1;
				}
				else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
				{
					int pp = elements[i].nodes_id[j];
					elements[i].nodes_id[j] = new_count[pp];
				}
				else
				{
					cout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for pure nonlocal fracture propagation!" << endl;	
					hout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for pure nonlocal fracture propagation!" << endl;
					return 0;
				}
			}
			if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
		}

		//--------------------------------------------
		//clear the relative elements of nodes
		for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_eles.clear();
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Read 3D mesh data from the file
int Mesher::Read_mesh(const string &output_file_name)
{
	ifstream idata(output_file_name.c_str());
	string s;
	for(int i=0; i<3; i++)	getline(idata, s);
	idata >> s;
	//---------------------------------------------------------------------------
	//Read nodes information
	vector<Node> temp_nod;
	if(s=="$Nodes")
	{
		int nodn;
		idata >> nodn;
		for(int i=0; i<nodn; i++)
		{
			int nnum;
			double x, y, z;
			idata >> nnum >> x >> y >> z;
			Node nd(x,y,z);
			nd.type = -1; //The type of node is unknown (0: internal node, 1: surface node, 2: boundary node; 3 corner node, -1: unknown)
			nd.DG_node=false;
			temp_nod.push_back(nd);
		}
		if(nodn!=(int)temp_nod.size()) 
		{
			cout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			hout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			return 0;
		}
	}
	else
	{
		cout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		hout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		return 0;
	}

	idata >> s;
	if(s!="$EndNodes") 
	{
		cout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		hout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		return 0;
	}

	//new number of nodes in vector (delete useless nodes which were used for creating geometry)
	vector<int> num_nod(temp_nod.size(), -1);
	//---------------------------------------------------------------------------
	//Read elements information
	idata >> s;
	if(s=="$Elements")
	{
		int elen;
		idata >> elen;
		for(int i=0; i<elen; i++)
		{
			int elnum, etype;
			string temps;
			idata >> elnum >> etype;

			if(etype!=5) getline(idata, temps);
			else
			{
				int temd;  
				for(int j=0; j<3; j++) idata >> temd;

				Element ele_temp;
				for(int j=0; j<8; j++)
				{
					idata >> temd;
					temd--;
					if(num_nod[temd]==-1)
					{
						nodes.push_back(temp_nod[temd]);
						num_nod[temd] = (int)nodes.size()-1;
					}
					ele_temp.nodes_id.push_back(num_nod[temd]);
				}
				ele_temp.DG_elem = false;
				ele_temp.type = 381;
				ele_temp.mat = 0;
				ele_temp.flag = 1;
				elements.push_back(ele_temp);
			}
		}
	}

	idata >> s;
	if(s!="$EndElements") 
	{
		cout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		hout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Read 2D mesh data from the file
int Mesher::Read_mesh_2DTet(const string &output_file_name)
{
	ifstream idata(output_file_name.c_str());
	string s;
	for(int i=0; i<3; i++)	getline(idata, s);
	idata >> s;
	//---------------------------------------------------------------------------
	//Read nodes information
	vector<Node> temp_nod;
	if(s=="$Nodes")
	{
		int nodn;
		idata >> nodn;
		for(int i=0; i<nodn; i++)
		{
			int nnum;
			double x, y, z;
			idata >> nnum >> x >> y >> z;
			Node nd(x,y,z);
			nd.type = -1;				//The type of node is unknown (0: internal node, 1: surface node, 2: boundary node; 3 corner node, -1: unknown)
			nd.DG_node=false;
			temp_nod.push_back(nd);
		}
		if(nodn!=(int)temp_nod.size()) 
		{
			cout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			hout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			return 0;
		}
	}
	else
	{
		cout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		hout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		return 0;
	}

	idata >> s;
	if(s!="$EndNodes") 
	{
		cout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		hout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		return 0;
	}

	//new number of nodes in vector (delete useless nodes which were used for creating geometry)
	vector<int> num_nod(temp_nod.size(), -1);
	//---------------------------------------------------------------------------
	//Read elements information
	idata >> s;
	if(s=="$Elements")
	{
		int elen;
		idata >> elen;
		for(int i=0; i<elen; i++)
		{
			int elnum, etype;
			string temps;
			idata >> elnum >> etype;

			if(etype!=3) getline(idata, temps);
			else
			{
				int temd;  
				for(int j=0; j<3; j++) idata >> temd;

				Element ele_temp;
				for(int j=0; j<4; j++)
				{
					idata >> temd;
					temd--;
					if(num_nod[temd]==-1)
					{
						nodes.push_back(temp_nod[temd]);
						num_nod[temd] = (int)nodes.size()-1;
					}
					ele_temp.nodes_id.push_back(num_nod[temd]);
				}
				ele_temp.DG_elem = false;
				ele_temp.type = 241;
				ele_temp.mat = 0;
				ele_temp.flag = 1;
				elements.push_back(ele_temp);
			}
		}
	}
	else
	{
		cout << "Error: the read string is " << s << ", not $Elements." << endl; 
		hout << "Error: the read string is " << s << ", not $Elements." << endl; 
		return 0;
	}

	idata >> s;
	if(s!="$EndElements") 
	{
		cout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		hout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Generate brick background grids
int Mesher::Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//read element size
	dx = grid_size.delta_x;
	dy = grid_size.delta_y;
	dz = grid_size.delta_z;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��������(�ڵ�͵�Ԫ)

	//���ɸ�������ĵ���
	int x_sec_num = int(geom_rve.len_x/dx+Zero) +1;
	int y_sec_num = int(geom_rve.wid_y/dy+Zero)+1;
	int z_sec_num = int(geom_rve.hei_z/dz+Zero) +1;

	//΢�����ʷֲ�������֤�����ĳߴ粻��
	dx = geom_rve.len_x/(x_sec_num-1);
	dy = geom_rve.wid_y/(y_sec_num-1);
	dz = geom_rve.hei_z/(z_sec_num-1);

	for( int i=0; i<z_sec_num; i++ )
	{
		double z = geom_rve.origin.z + i * dz ;
		for( int j=0; j<y_sec_num; j++ )
		{
			double y = geom_rve.origin.y + j * dy ;
			for( int k=0; k<x_sec_num; k++ )
			{
				double x = geom_rve.origin.x + k * dx ;

				Node nd(x,y,z);
				nd.type = Deter_node_type(i, j, k, z_sec_num-1, y_sec_num-1, x_sec_num-1);	//��ע�ڵ��λ��(�ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ�)
				nd.DG_node=false;
				nodes.push_back(nd);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���������嵥Ԫ

	//��������������
	for( int i=0; i<z_sec_num-1; i++ )
	{
		for( int j=0; j<y_sec_num-1; j++ )
		{
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//������İ˸�����
				Element ele_temp;
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
				ele_temp.DG_elem = false;

				elements.push_back(ele_temp);
			}
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//Generate brick grids for static fracture problems
int Mesher::Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve, const string &disc, const struct Weight_func &weight_func)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//read element size
	dx = grid_size.delta_x;
	dy = grid_size.delta_y;
	dz = grid_size.delta_z;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��������(�ڵ�͵�Ԫ)

	//���ɸ�������ĵ���
	int x_sec_num = int(geom_rve.len_x/dx+Zero) +1;
	int y_sec_num = int(geom_rve.wid_y/dy+Zero)+1;
	int z_sec_num = int(geom_rve.hei_z/dz+Zero) +1;

	//΢�����ʷֲ�������֤�����ĳߴ粻��
	dx = geom_rve.len_x/(x_sec_num-1);
	dy = geom_rve.wid_y/(y_sec_num-1);
	dz = geom_rve.hei_z/(z_sec_num-1);

	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		for( int i=0; i<z_sec_num; i++ )
		{
			double z = geom_rve.origin.z + i * dz ;
			for( int j=0; j<y_sec_num; j++ )
			{
				double y = geom_rve.origin.y + j * dy ;
				for( int k=0; k<x_sec_num; k++ )
				{
					double x = geom_rve.origin.x + k * dx ;

					Node nd(x,y,z);
					nd.type = Deter_node_type(i, j, k, z_sec_num-1, y_sec_num-1, x_sec_num-1);	//��ע�ڵ��λ��(�ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ�)
					nd.DG_node=false;
					nodes.push_back(nd);
				}
			}
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//���������嵥Ԫ

		//��������������
		for( int i=0; i<z_sec_num-1; i++ )
		{
			for( int j=0; j<y_sec_num-1; j++ )
			{
				for( int k=0; k<x_sec_num-1; k++ )
				{
					//������İ˸�����
					Element ele_temp;
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
					ele_temp.DG_elem = false;

					elements.push_back(ele_temp);
				}
			}
		}
	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		vector<Node> nodesFEM;
		vector<int> new_count;
		WeightFunc Wei_Fun;  //the class of weighting function
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//���ɽڵ�
		for( int i=0; i<z_sec_num; i++ )
		{
			double z = geom_rve.origin.z + i * dz ;
			for( int j=0; j<y_sec_num; j++ )
			{
				double y = geom_rve.origin.y + j * dy ;
				for( int k=0; k<x_sec_num; k++ )
				{
					double x = geom_rve.origin.x + k * dx ;
					Node nd(x,y,z);
					Point_3D poi_tem(nd.x, nd.y, nd.z);
					nd.type = Deter_node_type(i, j, k, z_sec_num-1, y_sec_num-1, x_sec_num-1);	//��ע�ڵ��λ��(�ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ�)
					if(fabs(Wei_Fun.Value_weight_function(poi_tem, weight_func)-1.0)<=Zero)
					{
						nd.DG_node = true;
						new_count.push_back(-1);
					}
					else
					{
						nd.DG_node = false;
						new_count.push_back((int)nodes.size());
						nodes.push_back(nd);
					}
					nodesFEM.push_back(nd);
				}
			}
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//���������嵥Ԫ

		//��������������
		for( int i=0; i<z_sec_num-1; i++ )
		{
			for( int j=0; j<y_sec_num-1; j++ )
			{
				for( int k=0; k<x_sec_num-1; k++ )
				{
					//������İ˸�����
					Element ele_temp;
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
					ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);

					ele_temp.DG_elem=false;
					int counter=0;
					for(int m=0; m<8; m++)
					{
						if(nodesFEM[ele_temp.nodes_id[m]].DG_node&&new_count[ele_temp.nodes_id[m]]==-1)
						{
							counter++;
							nodes.push_back(nodesFEM[ele_temp.nodes_id[m]]);
							ele_temp.nodes_id[m] = (int)nodes.size()-1;
						}
						else if(!nodesFEM[ele_temp.nodes_id[m]].DG_node&&new_count[ele_temp.nodes_id[m]]!=-1)
						{
							int pp = ele_temp.nodes_id[m];
							ele_temp.nodes_id[m] = new_count[pp];
						}
						else
						{
							cout << "Error: the specified situation is not defined in DGFEM part of \"Brick_background_mesh\" fucntion!" << endl;	
							hout << "Error: the specified situation is not defined in DGFEM part of \"Brick_background_mesh\" fucntion!" << endl;
							return 0;
						}
					}
					if(counter==8)	ele_temp.DG_elem=true;

					elements.push_back(ele_temp);
				}
			}
		}
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	return 1;

}
//---------------------------------------------------------------------------
//Generate brick grids for damage and fracture propagation
int Mesher::Brick_background_mesh(const struct Grid_size &grid_size, const struct Geom_RVE &geom_rve, const string &disc, const vector<bool> &full_dam_eles)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Brick_background_mesh(grid_size, geom_rve)== 0) return 0;

	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Brick_background_mesh(grid_size, geom_rve)== 0) return 0;
		
		//-----------------------------------------------------------------------------------------------------------------------------------------
		if(full_dam_eles.size()==0) return 1;
		
		//-----------------------------------------------------------------------------------------------------------------------------------------
		Deter_nodes_relative_eles();

		//-----------------------------------------------------------------------------------------------------------------------------------------
		vector<Node> nodesFEM(nodes);
		nodes.clear();
		vector<int> new_count;
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Define the DG nodes
		for(int i=0; i<(int)nodesFEM.size(); i++)
		{
			//-----------------------------------------------------------------------
			//Damaged node, at least one of its relative elements is damaged.
			//for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
			//{
			//	if(full_dam_eles[nodesFEM[i].relative_eles[j]]) 
			//	{
			//		nodesFEM[i].DG_node = true;
			//		break;
			//	}
			//}

			//-----------------------------------------------------------------------
			//Full damaged node, Its relative elements are all damaged. 
			nodesFEM[i].DG_node = true;
			for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
			{
				if(!full_dam_eles[nodesFEM[i].relative_eles[j]]) 
				{
					nodesFEM[i].DG_node = false;
					break;
				}
			}			
			
			if(nodesFEM[i].DG_node) new_count.push_back(-1);
			else
			{
				new_count.push_back((int)nodes.size());
				nodes.push_back(nodesFEM[i]);
			}
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].DG_elem=false;
			int counter=0;
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
				{
					counter++;
					nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
					elements[i].nodes_id[j] = (int)nodes.size()-1;
				}
				else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
				{
					int pp = elements[i].nodes_id[j];
					elements[i].nodes_id[j] = new_count[pp];
				}
				else
				{
					cout << "Error: the specified situation is not defined in DGFEM part of \"Brick_background_mesh\" fucntion for damage and fracture propagation!" << endl;	
					hout << "Error: the specified situation is not defined in DGFEM part of \"Brick_background_mesh\" fucntion for damage and fracture propagation!" << endl;
					return 0;
				}
			}
			if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//clear the relative elements of nodes
		for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_eles.clear();
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	return 1;
}
//---------------------------------------------------------------------------
//��������Ԫ��������
int Mesher::Element_material_property(const struct Ele_prop &ele_prop)
{	
	if(ele_prop.type=="Pure")	//pure material
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;		//��Ӧ��Ԫ����״����(��������xyz�� x ά�ȣ�y��Ԫ�Ľڵ������z��Ԫ���κ����ݴ�)��
														//���磬121: һά���ڵ�(�߶�)�����κ�����231: ��ά���ڵ�(������)�����κ�����241: ��ά�Ľڵ�(�ı���)�����κ�����
														//341: ��ά�Ľڵ�(������)�����κ�����361: ��ά���ڵ�(������)�����κ�����381����ά�˽ڵ�(������)�����κ���
			elements[i].flag = 0;			//��Ӧ��Ԫ�仯��ǣ�0 ��ʼ��̬����׼�����壩��1�仯����̬���������壬���Ǳ�׼�����壩
			elements[i].mat = 0;			//the type of material properties for damage part
		}
	}
	else if(ele_prop.type=="Fiber") //unidirectional fiber
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;	
			elements[i].flag = 0;			//��Ӧ��Ԫ�仯��ǣ�0 ��ʼ��̬����׼�����壩��1�仯����̬���������壬���Ǳ�׼�����壩
			Point_3D centp(0,0,0);
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				centp.x += nodes[elements[i].nodes_id[j]].x;
				centp.y += nodes[elements[i].nodes_id[j]].y;
				centp.z += nodes[elements[i].nodes_id[j]].z;
			}
			centp = centp/elements[i].nodes_id.size();
			if((centp.x-0.5)*(centp.x-0.5)+(centp.y-0.5)*(centp.y-0.5)<ele_prop.radius*ele_prop.radius) elements[i].mat = 1;
			else elements[i].mat = 0;
		}
	}
	else if(ele_prop.type=="Zone")
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;	
			elements[i].flag = 0;			//��Ӧ��Ԫ�仯��ǣ�0 ��ʼ��̬����׼�����壩��1�仯����̬���������壬���Ǳ�׼�����壩
			Point_3D centp(0,0,0);
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				centp.x += nodes[elements[i].nodes_id[j]].x;
				centp.y += nodes[elements[i].nodes_id[j]].y;
				centp.z += nodes[elements[i].nodes_id[j]].z;
			}
			centp = centp/elements[i].nodes_id.size();
			if(centp.x>ele_prop.zone_xmin&&centp.x<ele_prop.zone_xmax&&
				centp.y>ele_prop.zone_ymin&&centp.y<ele_prop.zone_ymax&&
				centp.z>ele_prop.zone_zmin&&centp.z<ele_prop.zone_zmax)
			{
				elements[i].mat = 1;
			}
			else
			{
				elements[i].mat = 0;
			}
		}
	}
	else
	{ 
		cout <<"Error: the specified type of material "<< ele_prop.type <<" is not defined yet!" << endl;	
		hout <<"Error: the specified type of material "<< ele_prop.type <<" is not defined yet!" << endl;	
		return 0;
	}
	return 1;
}
//---------------------------------------------------------------------------
//���Tecplot���ӻ���������
int Mesher::Export_mesh_data_tecplot(string output_file_name)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Tecplot_Meshes" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for (int i=0; i < (int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close();
	return 1;
}
//---------------------------------------------------------------------------
//���ݸ�����i j k�Լ����i_max j_max k_max�����������ڵ��λ�ã��ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ㣩
//����ֵ��0�ڵ㣬1�߽���㣬2�߽��ߵ㣬3�ǵ�
//i--z���ꣻ j--y���ꣻk--x����
int Mesher::Deter_node_type(const int i, const int j, const int k, const int i_max, const int j_max, const int k_max)const
{

	int type=0;
	if( i == 0 )	type++;
	else if( i == i_max ) type++;

	if( j == 0 )	type++;
	else if( j == j_max ) type++;

	if( k == 0 )	type++;
	else if( k == k_max ) type++;

	return type;
}
//---------------------------------------------------------------------------
//���ɸ�������Ԥ���ֲ�ģ�͵�Чģ��
int Mesher::Generate_grids_for_effective_stiffness(const struct Grid_size &nonloc_gsize, const double &decayR, double &grid_vol)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������ӳߴ�(������͸�)
	dx=nonloc_gsize.delta_x;
	dy=nonloc_gsize.delta_y;
	dz=nonloc_gsize.delta_z;

	//���ӵ����
	grid_vol = dx*dy*dz;

	//���ɸ�������ĵ���
	int x_sec_num = 2*(int(decayR/dx+Zero)+1);  //��һ��СֵΪ�˱���ȡ���Ĵ���
	int y_sec_num = 2*(int(decayR/dy+Zero)+1);
	int z_sec_num = 2*(int(decayR/dy+Zero)+1);

	Point_3D origin(0,0,0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɽڵ�
	for( int i=0; i<z_sec_num; i++ )
	{
		double z = origin.z + i * dz ;
		for( int j=0; j<y_sec_num; j++ )
		{
			double y = origin.y + j * dy ;
			for( int k=0; k<x_sec_num; k++ )
			{
				double x = origin.x + k * dx ;

				Node nd(x,y,z);
				nodes.push_back(nd);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���������嵥Ԫ

	//��������������
	for( int i=0; i<z_sec_num-1; i++ )
		for( int j=0; j<y_sec_num-1; j++ )
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//������İ˸�����
				Element ele_temp;
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);

				elements.push_back(ele_temp);
			}

	return 1;
}
//---------------------------------------------------------------------------
//���ɸ�������Ԥ���ֲ�ģ�͵�Чģ��(2D), z��ڵ�ֵ�����0
int Mesher::Generate_grids_for_effective_stiffness_2D(const struct Grid_size &nonloc_gsize, const double &decayR, double &grid_vol)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������ӳߴ�(���Ϳ�)
	dx=nonloc_gsize.delta_x;
	dy=nonloc_gsize.delta_y;

	//���ӵ����
	grid_vol = dx*dy;

	//���ɸ�������ĵ���
	int x_sec_num = 2*(int(decayR/dx+Zero)+1);  //��һ��СֵΪ�˱���ȡ���Ĵ���
	int y_sec_num = 2*(int(decayR/dy+Zero)+1);

	Point_3D origin(0,0,0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɽڵ�
	double z = origin.z;
	for( int i=0; i<y_sec_num; i++ )
	{
		double y = origin.y + i * dy ;
		for( int j=0; j<x_sec_num; j++ )
		{
			double x = origin.x + j * dx ;

			Node nd(x,y,z);
			nodes.push_back(nd);
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�����ı��ε�Ԫ
	for( int i=0; i<y_sec_num-1; i++ )
		for( int j=0; j<x_sec_num-1; j++ )
		{
			Element ele_temp;
			ele_temp.nodes_id.push_back(i * x_sec_num + j);
			ele_temp.nodes_id.push_back(i * x_sec_num + j + 1);
			ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num + j + 1);
			ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num + j);

			elements.push_back(ele_temp);
		}

	return 1;
}
//---------------------------------------------------------------------------
//����������������
int Mesher::Add_crack_to_element(const struct Crack &Cracks, const struct Geom_RVE &geom_rve)
{
	//loop by the number of cracks
	for(int i=0; i<Cracks.num; i++)
	{
		//---------------------------------------------------------------------------
		//reads crack parameteres
		struct Crack_Paras crapa;		//structure that contains crack parameters	
		if(Single_crack_parameters(i, Cracks, geom_rve, crapa)==0) return 0;

		//---------------------------------------------------------------------------
		//define relative elements to nodes
		if(Deter_nodes_relative_eles()==0 ) return 0;

		//---------------------------------------------------------------------------
		//�������Ƶ�������
		if(Adjust_mesh_with_crack(crapa)==0) return 0; 

		//---------------------------------------------------------------------------
		//�������ƺ�����ڵ����ص�Ԫ��Ϣ
		for(int j=0; j<(int)nodes.size(); j++)	nodes[j].relative_eles.clear(); //�������ڵ����ص�Ԫ��Ϣ
	}

	return 1;
}
//---------------------------------------------------------------------------
//the parameters of single crack
int Mesher::Single_crack_parameters(const int &nc, const struct Crack &Cracks, const struct Geom_RVE &geom_rve, struct Crack_Paras &crapa)
{
	crapa.num = nc;
	crapa.n0 = Cracks.n0[nc];		crapa.ch0 = Cracks.ch0[nc];
	crapa.n1 = Cracks.n1[nc];		crapa.ch1 = Cracks.ch1[nc];		crapa.d10 = Cracks.d10[nc];
	crapa.n2 = Cracks.n2[nc];		crapa.ch2 = Cracks.ch2[nc];		crapa.d20 = Cracks.d20[nc];		crapa.d21 = Cracks.d21[nc];
	crapa.n3 = Cracks.n3[nc];		crapa.d30 = Cracks.d30[nc];		crapa.d31 = Cracks.d31[nc];

	if((crapa.ch1=='x'&&(crapa.d10<=geom_rve.origin.x||crapa.d10>=geom_rve.origin.x+geom_rve.len_x))||
		(crapa.ch1=='y'&&(crapa.d10<=geom_rve.origin.y||crapa.d10>=geom_rve.origin.y+geom_rve.wid_y))||
		(crapa.ch1=='z'&&(crapa.d10<=geom_rve.origin.z||crapa.d10>=geom_rve.origin.z+geom_rve.hei_z))) 
	{

		cout << "Error: The parameter in the line 1 of  the crack " << nc << " is defined incorrectly!" << endl;
		hout << "Error: The parameter in the line 1 of  the crack " << nc << " is defined incorrectly!" << endl;
		return 0;
	}

	if((crapa.ch2=='x'&&(crapa.d21<=geom_rve.origin.x||crapa.d20>=geom_rve.origin.x+geom_rve.len_x))||
		(crapa.ch2=='y'&&(crapa.d21<=geom_rve.origin.y||crapa.d20>=geom_rve.origin.y+geom_rve.wid_y))||
		(crapa.ch2=='z'&&(crapa.d21<=geom_rve.origin.z||crapa.d20>=geom_rve.origin.z+geom_rve.hei_z))) 
	{

		cout << "Error: The parameter in the line 2 of  the crack " << nc << " is defined incorrectly!" << endl;
		hout << "Error: The parameter in the line 2 of  the crack " << nc << " is defined incorrectly!" << endl;
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//ȷ�����нڵ����ص�Ԫ��Ϣ
int Mesher::Deter_nodes_relative_eles()
{
	for( int i=0; i<(int)nodes.size(); i++ )
	{
		nodes[i].relative_eles.clear();
	}
	for( int i=0; i<(int)elements.size(); i++ )
	{
		int node_size = int(elements[i].nodes_id.size());
		for( int j=0; j<node_size; j++ )
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}
	return 1;
}
//---------------------------------------------------------------------------
//�������Ƶ�������
int Mesher::Adjust_mesh_with_crack(const struct Crack_Paras &crapa)
{
	int NS = (int)nodes.size();	//��Ϊ�����¹�����Ҫ����½ڵ�
	for(int i=0; i<NS; i++)
	{
		if(nodes[i].DG_node) continue;
		if((crapa.ch1=='x'&&fabs(nodes[i].x-crapa.d10)<Zero)||
			(crapa.ch1=='y'&&fabs(nodes[i].y-crapa.d10)<Zero)||
			(crapa.ch1=='z'&&fabs(nodes[i].z-crapa.d10)<Zero))		//�õ���ƽ����
		{
			if((crapa.ch2=='x'&&nodes[i].x>crapa.d20-Zero&&nodes[i].x<crapa.d21+Zero)||
				(crapa.ch2=='y'&&nodes[i].y>crapa.d20-Zero&&nodes[i].y<crapa.d21+Zero)||
				(crapa.ch2=='z'&&nodes[i].z>crapa.d20-Zero&&nodes[i].z<crapa.d21+Zero))	//�õ����߶���
			{
				double point_distance=0;
				if(crapa.ch2=='x') point_distance = fabs(nodes[i].x-crapa.d20);
				else if(crapa.ch2=='y') point_distance = fabs(nodes[i].y-crapa.d20);
				else if(crapa.ch2=='z') point_distance = fabs(nodes[i].z-crapa.d20);

				double crack_line_length = fabs(crapa.d21-crapa.d20);

				//����õ���ƶ�����
				double crack_t = crapa.d30;			//λ�õı���ֵ[0,1]
				double crack_width = 0;
				if(crapa.ch1=='x') crack_width = crapa.d31*dx;
				else if(crapa.ch1=='y') crack_width = crapa.d31*dy;
				else if(crapa.ch1=='z') crack_width = crapa.d31*dz;

				double distance_traveled = 0.0;
				double compare_val = point_distance-crack_t*crack_line_length;
				if(fabs(compare_val)<=Zero) distance_traveled = 0.5*crack_width;
				else if(compare_val<-Zero) distance_traveled = (point_distance/(crack_t*crack_line_length))*0.5*crack_width;
				else if(compare_val>Zero) distance_traveled = ((crack_line_length-point_distance)/((1-crack_t)*crack_line_length))*0.5*crack_width;

				//����ƶ��������0���Ͳ����ƶ��ˣ�����������ܳ����������˵�
				if(fabs(distance_traveled)<Zero)  continue;

				//�����½ڵ�
				Node new_node(nodes[i].x, nodes[i].y, nodes[i].z);
				new_node.type = nodes[i].type;
				//�޸Ľڵ�����
				if(crapa.ch1=='x')
				{
					new_node.x = nodes[i].x + distance_traveled;
					nodes[i].x = nodes[i].x - distance_traveled;
				}
				else if(crapa.ch1=='y')
				{
					new_node.y = nodes[i].y + distance_traveled;
					nodes[i].y = nodes[i].y - distance_traveled;
				}
				else if(crapa.ch1=='z')
				{
					new_node.z = nodes[i].z + distance_traveled;
					nodes[i].z = nodes[i].z - distance_traveled;
				}

				//�����½ڵ㵽�ڵ�����
				nodes.push_back(new_node);
				int new_node_num = (int)nodes.size()-1;
				//�޸���ص�Ԫ��Ϣ
				for(int j=0; j<(int)nodes[i].relative_eles.size(); j++)
				{
					double value =0;
					int ntrev = nodes[i].relative_eles[j];
					if(elements[ntrev].flag==0) elements[ntrev].flag = 1; //��Ԫ����״�б仯
					int etnis = (int)elements[ntrev].nodes_id.size();
					int nik = 0;		//���ڼ�¼�ڵ��ڵ�Ԫ�еľֲ����
					int nik_count = 0;		//���ڼ�¼����ص�Ԫ���Ƿ������һ���˽ڵ�
					for(int k=0; k<etnis; k++)
					{
						if(elements[ntrev].nodes_id[k] == i) { nik = k; nik_count++; continue; }		//��¼�ڵ��ڵ�Ԫ�еľֲ���Ų���¼����ص�Ԫ���Ƿ������һ���˽ڵ�
						if(crapa.ch1=='x')	value += nodes[elements[ntrev].nodes_id[k]].x;
						else if(crapa.ch1=='y') value += nodes[elements[ntrev].nodes_id[k]].y;
						else if(crapa.ch1=='z') value += nodes[elements[ntrev].nodes_id[k]].z;		
					}
					if(nik_count!=1)
					{
						hout <<"��" << ntrev<< "�ŵ�Ԫ�а���"<<nik_count<<"��"<< i <<"�Žڵ㡣" << endl;
						return 0;
					}
					value = value/(etnis-1);
					//�ж����ĵ��Ƿ�͵�i�ŵ���࣬������˵�����²���Ľڵ�ͬ��
					if(value>crapa.d10)
					{
						elements[ntrev].nodes_id[nik] = new_node_num;
						nodes[new_node_num].relative_eles.push_back(ntrev);  //���뵽�½ڵ����ص�Ԫ������
						nodes[i].relative_eles.erase(nodes[i].relative_eles.begin()+j);  //�ڵ�i�ŵ����ص�Ԫ��ɾ���˽ڵ�
						j--;
					}
				}
			}
		}
	}
	return 1;
}
//---------------------------------------------------------------------------
//Defines the neigbours of element and node (New Version with weight_func determination)
int Mesher::Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks, const struct Weight_func &weight_func)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Indicates itself in relative elements
	const int ES = (int)elements.size();
	if(mod=="Local")
	{
		for(int i=0; i<ES; i++) 
			elements[i].relative_eles.push_back(i);
	}
	else if(mod=="Hybrid")
	{
		if(weight_func.num==1&&weight_func.shape[0]=="Null")
		{
			for(int i=0; i<ES; i++) 
				elements[i].relative_eles.push_back(i);
		}
		else
		{
			//Non-local element neighbours
			vector<Point_3D> centrele(elements.size());
			for(int i=0; i<ES; i++)
			{
				Point_3D poi_tem(0,0,0);
				int ensize = elements[i].nodes_id.size();
				for(int j=0; j<ensize; j++)
				{
					poi_tem.x += nodes[elements[i].nodes_id[j]].x;
					poi_tem.y += nodes[elements[i].nodes_id[j]].y;
					poi_tem.z += nodes[elements[i].nodes_id[j]].z;
				}
				centrele[i] = poi_tem/ensize;
			}

			for(int i=0; i<ES; i++)   
				for(int j=0; j<ES; j++)
				{
					bool mark = false;
					double twodist = 2*dist;
					if(fabs(centrele[j].x-centrele[i].x)<=dist&&fabs(centrele[j].y-centrele[i].y)<=dist&&fabs(centrele[j].z-centrele[i].z)<=dist) mark = true;
					else if(fabs(centrele[j].x-centrele[i].x)<=twodist&&fabs(centrele[j].y-centrele[i].y)<=twodist&&fabs(centrele[j].z-centrele[i].z)<=twodist)
					{
						int elsize = elements[i].nodes_id.size();
						int ersize = elements[j].nodes_id.size();
						for(int k=0; k<elsize; k++)
						{
							for(int m=0; m<ersize; m++)
							{
								if(fabs(nodes[elements[i].nodes_id[k]].x-nodes[elements[j].nodes_id[m]].x)<=dist&&
								   fabs(nodes[elements[i].nodes_id[k]].y-nodes[elements[j].nodes_id[m]].y)<=dist&&
								   fabs(nodes[elements[i].nodes_id[k]].z-nodes[elements[j].nodes_id[m]].z)<=dist) 
								{	
									mark = true;
									break;
								}
							}
							if(mark) break;
						}
					}
					if(mark)
					{
						//Evaluate the presence of non-local cracks according to crack parameters
						if(cracks.num==0) elements[i].relative_eles.push_back(j);
						else
						{
							for(int k=0; k<cracks.num; k++)
							{
								if(cracks.ch0[k]=='x'&&cracks.ch1[k]=='y')	Evaluate_relative_element_cracks(centrele[i].y, centrele[i].z, centrele[j].y, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if (cracks.ch0[k]=='x'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].y, centrele[j].z, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].z, centrele[j].x, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].x, centrele[j].z, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].y, centrele[j].x, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='y') Evaluate_relative_element_cracks(centrele[i].y, centrele[i].x, centrele[j].y, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							}
						}
					}
				}
		}
	}
	else
	{
		cout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		hout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//defines node relative neighbours
	for(int i=0; i<ES; i++)
	{
		int einsize = elements[i].nodes_id.size();
		for(int j=0; j<einsize; j++)
		{
			const int nodn = elements[i].nodes_id[j];
			const int eres = (int)elements[i].relative_eles.size();
			for(int k=0; k<eres; k++)
			{
				const int eiek = elements[i].relative_eles[k];
				int eknsize = elements[eiek].nodes_id.size();
				for(int m=0; m<eknsize; m++)
				{
					const int renod = elements[eiek].nodes_id[m];
					if(nodn!=renod)
					{
						//���ַ�����
						int left = 0;
						int right = (int)nodes[nodn].relative_nods.size()-1;
						while(right>=left)
						{
							int middle = (left + right)/2;
							if(nodes[nodn].relative_nods[middle] == renod) goto Node_Num_Same; //�ڵ�����ͬ�����
							else if(nodes[nodn].relative_nods[middle] > renod) right = middle - 1;
							else left = middle + 1;
						}
						nodes[nodn].relative_nods.insert(nodes[nodn].relative_nods.begin()+left, renod);
Node_Num_Same: ;					
					}
				}
			}
		}
	}
	return 1;
}
//---------------------------------------------------------------------------
//Find the neighbour elements of every element
int Mesher::Deter_relative_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks, const struct Weight_func &weight_func)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Indicates itself in relative elements
	const int ES = (int)elements.size();
	if(mod=="Local")
	{
		for(int i=0; i<ES; i++) 
			elements[i].relative_eles.push_back(i);
	}
	else if(mod=="Hybrid")
	{
		if(weight_func.num==1&&weight_func.shape[0]=="Null")
		{
			for(int i=0; i<ES; i++) 
				elements[i].relative_eles.push_back(i);
		}
		else
		{
			//Non-local element neighbours
			vector<Point_3D> centrele(elements.size());
			for(int i=0; i<ES; i++)
			{
				Point_3D poi_tem(0,0,0);
				int ensize = elements[i].nodes_id.size();
				for(int j=0; j<ensize; j++)
				{
					poi_tem.x += nodes[elements[i].nodes_id[j]].x;
					poi_tem.y += nodes[elements[i].nodes_id[j]].y;
					poi_tem.z += nodes[elements[i].nodes_id[j]].z;
				}
				centrele[i] = poi_tem/ensize;
			}

			for(int i=0; i<ES; i++)   
				for(int j=0; j<ES; j++)
				{
					bool mark = false;
					double twodist = 2*dist;
					if(fabs(centrele[j].x-centrele[i].x)<=dist&&fabs(centrele[j].y-centrele[i].y)<=dist&&fabs(centrele[j].z-centrele[i].z)<=dist) mark = true;
					else if(fabs(centrele[j].x-centrele[i].x)<=twodist&&fabs(centrele[j].y-centrele[i].y)<=twodist&&fabs(centrele[j].z-centrele[i].z)<=twodist)
					{
						int elsize = elements[i].nodes_id.size();
						int ersize = elements[j].nodes_id.size();
						for(int k=0; k<elsize; k++)
						{
							for(int m=0; m<ersize; m++)
							{
								if(fabs(nodes[elements[i].nodes_id[k]].x-nodes[elements[j].nodes_id[m]].x)<=dist&&
								   fabs(nodes[elements[i].nodes_id[k]].y-nodes[elements[j].nodes_id[m]].y)<=dist&&
								   fabs(nodes[elements[i].nodes_id[k]].z-nodes[elements[j].nodes_id[m]].z)<=dist) 
								{	
									mark = true;
									break;
								}
							}
							if(mark) break;
						}
					}
					if(mark)
					{
						//Evaluate the presence of non-local cracks according to crack parameters
						if(cracks.num==0) elements[i].relative_eles.push_back(j);
						else
						{
							for(int k=0; k<cracks.num; k++)
							{
								if(cracks.ch0[k]=='x'&&cracks.ch1[k]=='y')	Evaluate_relative_element_cracks(centrele[i].y, centrele[i].z, centrele[j].y, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if (cracks.ch0[k]=='x'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].y, centrele[j].z, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].z, centrele[j].x, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].x, centrele[j].z, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].y, centrele[j].x, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
								else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='y') Evaluate_relative_element_cracks(centrele[i].y, centrele[i].x, centrele[j].y, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							}
						}
					}
				}
		}
	}
	else
	{
		cout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		hout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		return 0;
	}
	return 1;
}
//---------------------------------------------------------------------------
//Defines the neigbours of element and node 
int Mesher::Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const string &mod, const struct Crack &cracks)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Indicates itself in relative elements
	const int ES = (int)elements.size();
	if(mod=="Local")
	{
		for(int i=0; i<ES; i++) 
			elements[i].relative_eles.push_back(i);
	}
	else if(mod=="Hybrid")
	{
		//Non-local element neighbours
		vector<Point_3D> centrele(elements.size());
		for(int i=0; i<ES; i++)
		{
			Point_3D poi_tem(0,0,0);
			int ensize = elements[i].nodes_id.size();
			for(int j=0; j<ensize; j++)
			{
				poi_tem.x += nodes[elements[i].nodes_id[j]].x;
				poi_tem.y += nodes[elements[i].nodes_id[j]].y;
				poi_tem.z += nodes[elements[i].nodes_id[j]].z;
			}
			centrele[i] = poi_tem/ensize;
		}

		for(int i=0; i<ES; i++)   
			for(int j=0; j<ES; j++)
			{
				bool mark = false;
				double twodist = 2*dist;
				if(fabs(centrele[j].x-centrele[i].x)<=dist&&fabs(centrele[j].y-centrele[i].y)<=dist&&fabs(centrele[j].z-centrele[i].z)<=dist) mark = true;
				else if(fabs(centrele[j].x-centrele[i].x)<=twodist&&fabs(centrele[j].y-centrele[i].y)<=twodist&&fabs(centrele[j].z-centrele[i].z)<=twodist)
				{
					int elsize = elements[i].nodes_id.size();
					int ersize = elements[j].nodes_id.size();
					for(int k=0; k<elsize; k++)
					{
						for(int m=0; m<ersize; m++)
						{
							if(fabs(nodes[elements[i].nodes_id[k]].x-nodes[elements[j].nodes_id[m]].x)<=dist&&
								fabs(nodes[elements[i].nodes_id[k]].y-nodes[elements[j].nodes_id[m]].y)<=dist&&
								fabs(nodes[elements[i].nodes_id[k]].z-nodes[elements[j].nodes_id[m]].z)<=dist) 
							{	
								mark = true;
								break;
							}
						}
						if(mark) break;
					}
				}
				if(mark)
				{
					//Evaluate the presence of non-local cracks according to crack parameters
					if(cracks.num==0) elements[i].relative_eles.push_back(j);
					else
					{
						for(int k=0; k<cracks.num; k++)
						{
							if(cracks.ch0[k]=='x'&&cracks.ch1[k]=='y')	Evaluate_relative_element_cracks(centrele[i].y, centrele[i].z, centrele[j].y, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							else if (cracks.ch0[k]=='x'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].y, centrele[j].z, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].z, centrele[j].x, centrele[j].z, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							else if(cracks.ch0[k]=='y'&&cracks.ch1[k]=='z') Evaluate_relative_element_cracks(centrele[i].z, centrele[i].x, centrele[j].z, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='x') Evaluate_relative_element_cracks(centrele[i].x, centrele[i].y, centrele[j].x, centrele[j].y, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
							else if(cracks.ch0[k]=='z'&&cracks.ch1[k]=='y') Evaluate_relative_element_cracks(centrele[i].y, centrele[i].x, centrele[j].y, centrele[j].x, cracks.d10[k], cracks.d20[k], cracks.d21[k], elements[i], j);
						}
					}
				}
			}
	}
	else
	{
		cout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		hout << "Error: the type of model (Model_Discretization) is wrong!" << endl;
		return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//defines node relative neighbours
	for(int i=0; i<ES; i++)
	{
		int einsize = elements[i].nodes_id.size();
		for(int j=0; j<einsize; j++)
		{
			const int nodn = elements[i].nodes_id[j];
			const int eres = (int)elements[i].relative_eles.size();
			for(int k=0; k<eres; k++)
			{
				const int eiek = elements[i].relative_eles[k];
				int eknsize = elements[eiek].nodes_id.size();
				for(int m=0; m<eknsize; m++)
				{
					const int renod = elements[eiek].nodes_id[m];
					if(nodn!=renod)
					{
						//���ַ�����
						int left = 0;
						int right = (int)nodes[nodn].relative_nods.size()-1;
						while(right>=left)
						{
							int middle = (left + right)/2;
							if(nodes[nodn].relative_nods[middle] == renod) goto Node_Num_Same; //�ڵ�����ͬ�����
							else if(nodes[nodn].relative_nods[middle] > renod) right = middle - 1;
							else left = middle + 1;
						}
						nodes[nodn].relative_nods.insert(nodes[nodn].relative_nods.begin()+left, renod);
Node_Num_Same: ;					
					}
				}
			}
		}
	}
	return 1;
}
//---------------------------------------------------------------------------
//Evaluate the relative elements under a crack condition
void Mesher::Evaluate_relative_element_cracks(const double &vi1, const double &vi2, const double &vj1, const double &vj2, const double &d10, const double &d20, const double &d21, Element &ele, const int &nj)const
{
	if((vi1-d10)*(vj1-d10)<=0)
	{
		double test=(vj2-vi2)/(vj1-vi1)*(d10-vi1)+vi2;
		if(test<=d20||test>=d21)	ele.relative_eles.push_back(nj);
	}
	else ele.relative_eles.push_back(nj);
}
//-----------------------------------------------------------------------------------------------
//Estimate the volume of brick element (deformed shape)
double Mesher::Calculate_brick_volume(const Node elenod[])const
{
	double sum = 0;
	int nod[6][4] = { {0,1,2,5}, {2,5,6,0}, {0,6,2,7}, {0,2,3,7}, {0,4,6,7}, {0,4,5,6} };  //������ָ������С������
	for(int tet=0; tet<6; tet++)
	{
		double a[4];
		for(int j=0; j<4; j++)
		{	
			a[j]	=	elenod[nod[tet][1]].x*(elenod[nod[tet][2]].y*elenod[nod[tet][3]].z-elenod[nod[tet][2]].z*elenod[nod[tet][3]].y)+
						elenod[nod[tet][1]].y*(elenod[nod[tet][2]].z*elenod[nod[tet][3]].x-elenod[nod[tet][2]].x*elenod[nod[tet][3]].z)+
						elenod[nod[tet][1]].z*(elenod[nod[tet][2]].x*elenod[nod[tet][3]].y-elenod[nod[tet][2]].y*elenod[nod[tet][3]].x);
			int temp = nod[tet][0];
			for(int k=0; k<3; k++) nod[tet][k] = nod[tet][k+1];
			nod[tet][3] = temp;
		}
		sum += a[0]-a[1]+a[2]-a[3];

	}
	sum = sum/6.0;

	return sum;
}
//-----------------------------------------------------------------------------------------------
//Estimate the volume of brick element (deformed shape)
double Mesher::Calculate_brick_volume(const double (*elenod)[3])const
{
	double sum = 0;
	int nod[6][4] = { {0,1,2,5}, {2,5,6,0}, {0,6,2,7}, {0,2,3,7}, {0,4,6,7}, {0,4,5,6} };  //������ָ������С������
	for(int tet=0; tet<6; tet++)
	{
		double a[4];
		for(int j=0; j<4; j++)
		{	
			a[j]	=	elenod[nod[tet][1]][0]*(elenod[nod[tet][2]][1]*elenod[nod[tet][3]][2]-elenod[nod[tet][2]][2]*elenod[nod[tet][3]][1])+
						elenod[nod[tet][1]][1]*(elenod[nod[tet][2]][2]*elenod[nod[tet][3]][0]-elenod[nod[tet][2]][0]*elenod[nod[tet][3]][2])+
						elenod[nod[tet][1]][2]*(elenod[nod[tet][2]][0]*elenod[nod[tet][3]][1]-elenod[nod[tet][2]][1]*elenod[nod[tet][3]][0]);
			int temp = nod[tet][0];
			for(int k=0; k<3; k++) nod[tet][k] = nod[tet][k+1];
			nod[tet][3] = temp;
		}
		sum += a[0]-a[1]+a[2]-a[3];

	}
	sum = sum/6.0;

	return sum;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh data(DGFEM) for 2D
int Mesher::Import_mesh_reconfiguration_2D(string output_file_name, const string disc, const vector<bool> &full_dam_eles, const struct Weight_func &weight_func)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

//		if(Read_mesh_3To2D(output_file_name, 0)== 0) return 0;
		if(Read_mesh_2DTet(output_file_name)== 0) return 0;

	}
	else if(disc=="DGFEM")
	{
		if(weight_func.num==1&&weight_func.shape[0]=="Null")
		{
			cout<< "    I'm using FEM elements and nodes."<<endl;
			hout<< "    I'm using FEM elements and nodes."<<endl;

			if(Read_mesh_2DTet(output_file_name)== 0) return 0;
		}
		else
		{
			cout<< "    I'm using DGFEM elements and nodes."<<endl;
			hout<< "    I'm using DGFEM elements and nodes."<<endl;

//			if(Read_mesh_3To2D(output_file_name, 0)== 0) return 0;
			if(Read_mesh_2DTet(output_file_name)== 0) return 0;

			//-----------------------------------------------------------------------------------------------------------------------------------------
			if(full_dam_eles.size()!=0)
			{
				//------------------------------------------
				Deter_nodes_relative_eles();

				//------------------------------------------
				vector<Node> nodesFEM(nodes);
				nodes.clear();
				vector<int> new_count;
				//------------------------------------------
				//Define the DG nodes
				for(int i=0; i<(int)nodesFEM.size(); i++)
				{
					//-----------------------------------------------------------------------
					//Damaged node, at least one of its relative elements is damaged.
					for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
					{
						if(full_dam_eles[nodesFEM[i].relative_eles[j]]) 
						{
							nodesFEM[i].DG_node = true;
							break;
						}
					}

					//-----------------------------------------------------------------------
					//Full damaged node, Its relative elements are all damaged. 
					//nodesFEM[i].DG_node = true;
					//for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
					//{
					//	if(!full_dam_eles[nodesFEM[i].relative_eles[j]]) 
					//	{
					//		nodesFEM[i].DG_node = false;
					//		break;
					//	}
					//}

					//-----------------------------------------------------------------------
					if(nodesFEM[i].DG_node) new_count.push_back(-1);
					else
					{
						nodes.push_back(nodesFEM[i]);
						new_count.push_back((int)nodes.size()-1);
					}
				}

				//------------------------------------------
				for(int i=0; i<(int)elements.size(); i++)
				{
					elements[i].flag = 1;
					elements[i].DG_elem=false;
					int counter=0;
					for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
					{
						if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
						{
							counter++;
							nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
							elements[i].nodes_id[j] = (int)nodes.size()-1;
						}
						else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
						{
							int pp = elements[i].nodes_id[j];
							elements[i].nodes_id[j] = new_count[pp];
						}
						else
						{
							cout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;	
							hout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;
							return 0;
						}
					}
					if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
				}

				//--------------------------------------------
				//clear the relative elements of nodes
				for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_eles.clear();
			}
		}
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot_2D("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh (DGFEM) by weighting function for 2D
int Mesher::Mesh_reconfiguration_weight_func_2D(string output_file_name, const string disc, const struct Weight_func &weight_func)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;
	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;

		//------------------------------------------
		vector<Node> nodesFEM(nodes);
		nodes.clear();
		vector<int> new_count;
		WeightFunc Wei_Fun;  //the class of weighting function
		//------------------------------------------
		//Define the DG nodes
		for(int i=0; i<(int)nodesFEM.size(); i++)
		{
			Point_3D poi_tem(nodesFEM[i].x, nodesFEM[i].y, nodesFEM[i].z);
			if(fabs(Wei_Fun.Value_weight_function(poi_tem, weight_func)-1.0)<=Zero)
			{
				nodesFEM[i].DG_node = true;
				new_count.push_back(-1);
			}
			else
			{
				nodesFEM[i].DG_node = false;
				new_count.push_back((int)nodes.size());
				nodes.push_back(nodesFEM[i]);
			}			
		}

		//------------------------------------------
		//Define the DG elements
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].flag = 1;
			elements[i].DG_elem=false;
			int counter=0;
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
				{
					counter++;
					nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
					elements[i].nodes_id[j] = (int)nodes.size()-1;
				}
				else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
				{
					int nod_num = elements[i].nodes_id[j];
					elements[i].nodes_id[j] = new_count[nod_num];
				}
				else
				{
					cout << "Error: the specified situation is not defined in DGFEM part of \"Mesh_reconfiguration_weight_func_2D\" fucntion for damage and fracture propagation!" << endl;	
					hout << "Error: the specified situation is not defined in DGFEM part of \"Mesh_reconfiguration_weight_func_2D\" fucntion for damage and fracture propagation!" << endl;
					return 0;
				}
			}
			if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
		}		
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot_2D("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh data(DGFEM) by weigted element for 2D
int Mesher::Mesh_reconfiguration_weighted_eles_2D(string output_file_name, const string disc, const vector<bool> &weighted_eles)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;

	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		if(weighted_eles.size()!=0)
		{
			//------------------------------------------
			Deter_nodes_relative_eles();

			//------------------------------------------
			vector<Node> nodesFEM(nodes);
			nodes.clear();
			vector<int> new_count;
			//------------------------------------------
			//Define the DG nodes
			for(int i=0; i<(int)nodesFEM.size(); i++)
			{
				//-----------------------------------------------------------------------
				//Full weighted node, Its relative elements are all weighted. 
				nodesFEM[i].DG_node = true;
				for(int j=0; j<(int)nodesFEM[i].relative_eles.size(); j++)
				{
					if(!weighted_eles[nodesFEM[i].relative_eles[j]]) 
					{
						nodesFEM[i].DG_node = false;
						break;
					}
				}

				if(nodesFEM[i].DG_node) new_count.push_back(-1);
				else
				{
					new_count.push_back((int)nodes.size());
					nodes.push_back(nodesFEM[i]);
				}
			}

			//------------------------------------------
			for(int i=0; i<(int)elements.size(); i++)
			{
				elements[i].flag = 1;
				elements[i].DG_elem=false;
				int counter=0;
				for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
				{
					if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
					{
						counter++;
						nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
						elements[i].nodes_id[j] = (int)nodes.size()-1;
					}
					else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
					{
						int pp = elements[i].nodes_id[j];
						elements[i].nodes_id[j] = new_count[pp];
					}
					else
					{
						cout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;	
						hout << "Error: the specified situation is not defined in DGFEM part of \"Import_mesh_reconfiguration\" fucntion for damage and fracture propagation!" << endl;
						return 0;
					}
				}
				if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
			}

			//--------------------------------------------
			//clear the relative elements of nodes
			for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_eles.clear();
		}
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot_2D("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Import mesh data from a file and reconfiguration mesh (DGFEM) by weighting function and delta legnth for 2D
int Mesher::Mesh_reconfiguration_weight_delta_2D(string output_file_name, const string disc, const struct Weight_func &weight_func, const struct Peri_para &peri_para)
{
	if(disc=="FEM")
	{
		cout<< "    I'm using FEM elements and nodes."<<endl;
		hout<< "    I'm using FEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;
	}
	else if(disc=="DGFEM")
	{
		cout<< "    I'm using DGFEM elements and nodes."<<endl;
		hout<< "    I'm using DGFEM elements and nodes."<<endl;

		if(Read_mesh_2DTet(output_file_name)== 0) return 0;

		//------------------------------------------
		vector<Node> nodesFEM(nodes);
		nodes.clear();
		vector<int> new_count;
		WeightFunc Wei_Fun;  //the class of weighting function
		
		//------------------------------------------
		//Define the weight values of nodes
		vector<int> node_weight(nodesFEM.size(), false);		
		for(int i=0; i<(int)nodesFEM.size(); i++)
		{
			Point_3D poi_tem(nodesFEM[i].x, nodesFEM[i].y, nodesFEM[i].z);
			if(fabs(Wei_Fun.Value_weight_function(poi_tem, weight_func)-1.0)<=Zero)
			{
				node_weight[i] = true;
			}
		}

		double dist = peri_para.horizon_R;
		//------------------------------------------
		//Define the DG nodes
		for(int i=0; i<(int)nodesFEM.size(); i++)
		{
			bool temp_key = true;
			if(node_weight[i])
			{
				for(int j=0; j<(int)nodesFEM.size(); j++)
				{
//					if(i!=j&&nodesFEM[i].distance_to(nodesFEM[j])<dist+Zero)
					if(i!=j&&
						fabs(nodesFEM[j].x-nodesFEM[i].x)<dist+Zero&&
						fabs(nodesFEM[j].y-nodesFEM[i].y)<dist+Zero&&
						fabs(nodesFEM[j].z-nodesFEM[i].z)<dist+Zero)
					{
						if(!node_weight[j])
						{	
							temp_key = false;
							break;
						}
					}
				}
			}
			else
			{
				temp_key =false;
			}

			if(temp_key)
			{
				nodesFEM[i].DG_node = true;
				new_count.push_back(-1);
			}
			else
			{
				nodesFEM[i].DG_node = false;
				new_count.push_back((int)nodes.size());
				nodes.push_back(nodesFEM[i]);
			}			
		}

		//------------------------------------------
		//Define the DG elements
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].flag = 1;
			elements[i].DG_elem=false;
			int counter=0;
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				if(nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]==-1)
				{
					counter++;
					nodes.push_back(nodesFEM[elements[i].nodes_id[j]]);
					elements[i].nodes_id[j] = (int)nodes.size()-1;
				}
				else if(!nodesFEM[elements[i].nodes_id[j]].DG_node&&new_count[elements[i].nodes_id[j]]!=-1)
				{
					int nod_num = elements[i].nodes_id[j];
					elements[i].nodes_id[j] = new_count[nod_num];
				}
				else
				{
					cout << "Error: the specified situation is not defined in DGFEM part of \"Mesh_reconfiguration_weight_func_2D\" fucntion for damage and fracture propagation!" << endl;	
					hout << "Error: the specified situation is not defined in DGFEM part of \"Mesh_reconfiguration_weight_func_2D\" fucntion for damage and fracture propagation!" << endl;
					return 0;
				}
			}
			if(counter==(int)elements[i].nodes_id.size())	elements[i].DG_elem=true;
		}		
	}
	
	cout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	hout << "    There is " << nodes.size() << " nodes in this simulation." << endl;
	cout << "    There is " << elements.size() << " elements in this simulation." << endl;
	hout << "    There is " << elements.size() << " elements in this simulation." << endl;

	//���Tecplot���ӻ���������
	if(Export_mesh_data_tecplot_2D("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Read mesh data from the file
int Mesher::Read_mesh_3To2D(const string &output_file_name, const double &Const_Z)
{
	ifstream idata(output_file_name.c_str());
	string s;
	for(int i=0; i<3; i++)	getline(idata, s);
	idata >> s;
	//---------------------------------------------------------------------------
	//Read nodes information
	vector<Node> temp_nod;
	if(s=="$Nodes")
	{
		int nodn;
		idata >> nodn;
		for(int i=0; i<nodn; i++)
		{
			int nnum;
			double x, y, z;
			idata >> nnum >> x >> y >> z;
			Node nd(x,y,z);
			nd.type = -1; //The type of node is unknown (0: internal node, 1: surface node, 2: boundary node; 3 corner node, -1: unknown)
			nd.DG_node=false;
			temp_nod.push_back(nd);
		}
		if(nodn!=(int)temp_nod.size()) 
		{
			cout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			hout << "Error: the number of nodes in data file: " << nodn << " is not same to the size of vector nodes: " << (int)nodes.size() << " !" << endl; 
			return 0;
		}
	}
	else
	{
		cout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		hout << "Error: the read string is " << s << ", not $Nodes." << endl; 
		return 0;
	}

	idata >> s;
	if(s!="$EndNodes") 
	{
		cout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		hout << "Error: the read string is " << s << ", not $EndNodes." << endl; 
		return 0;
	}

	//new number of nodes in vector (delete useless nodes which were used for creating geometry)
	vector<int> num_nod(temp_nod.size(), -1);
	//---------------------------------------------------------------------------
	//Read elements information
	idata >> s;
	if(s=="$Elements")
	{
		int elen;
		idata >> elen;
		for(int i=0; i<elen; i++)
		{
			int elnum, etype;
			string temps;
			idata >> elnum >> etype;

			if(etype!=5) getline(idata, temps);
			else
			{
				int temd;  
				for(int j=0; j<3; j++) idata >> temd;

				Element ele_temp;
				for(int j=0; j<8; j++)
				{
					idata >> temd;
					temd--;
					if(temp_nod[temd].z==Const_Z)
					{
						if(num_nod[temd]==-1)
						{
							nodes.push_back(temp_nod[temd]);
							num_nod[temd] = (int)nodes.size()-1;
						}
						ele_temp.nodes_id.push_back(num_nod[temd]);
					}
				}
				if(ele_temp.nodes_id.size()==4)
				{
					ele_temp.DG_elem = false;
					ele_temp.type = 241;
					ele_temp.mat = 0;
					ele_temp.flag = 1;
					elements.push_back(ele_temp);
				}
			}
		}
	}

	idata >> s;
	if(s!="$EndElements") 
	{
		cout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		hout << "Error: the read string is " << s << ", not $EndElements." << endl; 
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Export mesh data file for 2D grids shown in Tecplot
int Mesher::Export_mesh_data_tecplot_2D(string output_file_name)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Tecplot_Meshes" << endl;
	otec << "VARIABLES = X, Y" << endl;

	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for (int i=0; i < (int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y  << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Estimate the area of quadrilateral element
double Mesher::Calculate_quadri_area(const Node elenod[])const
{
	//calculate area (Wiki chinese)
	double sa=(elenod[0].x-elenod[1].x)*(elenod[0].x-elenod[1].x)
					+(elenod[0].y-elenod[1].y)*(elenod[0].y-elenod[1].y)
					+(elenod[0].z-elenod[1].z)*(elenod[0].z-elenod[1].z);

	double sb=(elenod[1].x-elenod[2].x)*(elenod[1].x-elenod[2].x)
					+(elenod[1].y-elenod[2].y)*(elenod[1].y-elenod[2].y)
					+(elenod[1].z-elenod[2].z)*(elenod[1].z-elenod[2].z);

	double sc=(elenod[2].x-elenod[3].x)*(elenod[2].x-elenod[3].x)
					+(elenod[2].y-elenod[3].y)*(elenod[2].y-elenod[3].y)
					+(elenod[2].z-elenod[3].z)*(elenod[2].z-elenod[3].z);
	
	double sd=(elenod[3].x-elenod[0].x)*(elenod[3].x-elenod[0].x)
					+(elenod[3].y-elenod[0].y)*(elenod[3].y-elenod[0].y)
					+(elenod[3].z-elenod[0].z)*(elenod[3].z-elenod[0].z);

	double se=(elenod[0].x-elenod[2].x)*(elenod[0].x-elenod[2].x)
					+(elenod[0].y-elenod[2].y)*(elenod[0].y-elenod[2].y)
					+(elenod[0].z-elenod[2].z)*(elenod[0].z-elenod[2].z);

	double sf=(elenod[1].x-elenod[3].x)*(elenod[1].x-elenod[3].x)
					+(elenod[1].y-elenod[3].y)*(elenod[1].y-elenod[3].y)
					+(elenod[1].z-elenod[3].z)*(elenod[1].z-elenod[3].z);
					
	double s_area = 0.25*sqrt(4*se*sf-(sb+sd-sa-sc)*(sb+sd-sa-sc));

	return s_area;
}
//===========================================================================
