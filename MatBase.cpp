//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatBase.cpp
//OBJECTIVE:	Generate material database
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#include "MatBase.h"

//---------------------------------------------------------------------------
//Generate the material base with damage information for 2D
int MatBase::Generate_matbase_2D(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para)
{
	MatPro Total_mat(stif_nonloc);
	if(Total_mat.Generate_elas_matrix_2D()==0) return 0; //���ɲ��ϵĸնȾ���

	//---------------------------------------------------------------------------
	MatPro known_mat;
	if(damages.d_crit<1.0)
	{
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
			{
				if(i==j&&(i==2||i==4||i==5))
				{
					known_mat.elas_matrix[i][j] = Total_mat.elas_matrix[i][j];
				}
				else
				{
					known_mat.elas_matrix[i][j] = Total_mat.elas_matrix[i][j]*(1.0-damages.d_crit);
				}
			}
		known_mat.Get_ele_para_by_ela_matrix();
	}
	else
	{
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				known_mat.elas_matrix[i][j] = 0.0;

		known_mat.E11 = known_mat.E22 = known_mat.E33 = 0.0;
		known_mat.Nu12 = known_mat.Nu13 = known_mat.Nu23 = 0.0;
		known_mat.G12 = known_mat.G13 = known_mat.G23 = 0.0;
		mats_vec.push_back(known_mat);
	}

	//---------------------------------------------------------------------------
	//�޸Ĳ���ϵ�������ݷǾֲ�������
	if(peri_para.horizon_R>0.0&&damages.d_crit<1.0) 
	{
		MatPro mat(stif_nonloc.type);
		//---------------------------------------------------------------------------
		//���ɸ�������Ԥ���ֲ�ģ�͵�Ч�ն���
		Mesher grid;
		double grid_vol;
		if(grid.Generate_grids_for_effective_stiffness_2D(nonloc_gsize, peri_para.horizon_R, grid_vol)==0) return 0;

		//ȡ��˹�㼰��Ȩϵ��
		Gauss gau;
		if(gau.Generate_gauss_2D(gaupoi_num)==0) return 0;
	
		//Ԥ���ֲ�ģ�͵�Ч�ն���
		MathMatrix Km(3, 3);
		for(int i=0; i<3; i++)
		{
			double a[3] = {0};
			a[i] = 1.0;
			if(Estimate_stiffness_long_range_interaction_2D(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, a, mat, peri_para.horizon_R, peri_para.intrinsic_L, 0)==0) return 0;
			Km.element[0][i] = mat.elas_matrix[0][0];
			Km.element[1][i] = mat.elas_matrix[1][1];
			Km.element[2][i] = mat.elas_matrix[0][1];
		}

		//���¼���
		double Cm[3] = {	known_mat.elas_matrix[0][0],	known_mat.elas_matrix[1][1],	known_mat.elas_matrix[0][1] };

		MathMatrix Dm(3,3);
		Dm = (Km.Transpose()*Km).Inverse()*Km.Transpose();
		for(int i=0; i<3; i++)
		{
			peri_para.acoe[i] = 0.0;
			for(int j=0; j<3; j++) peri_para.acoe[i] += Dm.element[i][j]*Cm[j];
		}

		if(Estimate_stiffness_long_range_interaction_2D(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, peri_para.acoe, mat, peri_para.horizon_R, peri_para.intrinsic_L, 1)==0) return 0;

		//---------------------------------------------------------------------------
		mat.elas_matrix[2][2] = stif_nonloc.E33;
		mat.elas_matrix[4][4] = stif_nonloc.G23;
		mat.elas_matrix[5][5] = stif_nonloc.G13;

		//---------------------------------------------------------------------------
		//���ݵ��Ծ��󣬷�����ϲ���
		mat.Get_ele_para_by_ela_matrix();

		mats_vec.push_back(mat); //�����������
	}

	//---------------------------------------------------------------------------
	//The total siffness matrix (for damage part)
	mats_vec.push_back(Total_mat);

	//---------------------------------------------------------------------------
	//���ݵ��Ծ��󣬷�����ϲ���
	Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 0);
	Print_stiffness_matrix("Material_stiffness_matrix.dat", Total_mat, 1);

	//---------------------------------------------------------------------------
	//������в��ϵĸնȾ��󼰶�Ӧ������ģ�������ɱ�
	Print_stiffness_matrix("Material_stiffness_matrix.dat", peri_para.acoe, 1);

	return 1;
}
//---------------------------------------------------------------------------
//Generate the material base
int MatBase::Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, struct Peri_para &peri_para)
{
	//---------------------------------------------------------------------------
	MatPro known_mat(stif_nonloc);
	if(known_mat.Generate_elas_matrix()==0) return 0; //���ɲ��ϵĸնȾ���

	//�޸Ĳ���ϵ�������ݷǾֲ�������
	if(peri_para.horizon_R>0.0) 
	{
		MatPro mat(stif_nonloc.type);
		//---------------------------------------------------------------------------
		//���ɸ�������Ԥ���ֲ�ģ�͵�Ч�ն���
		Mesher grid;
		double grid_vol;
		if(grid.Generate_grids_for_effective_stiffness(nonloc_gsize, peri_para.horizon_R, grid_vol)==0) return 0;

		//ȡ��˹�㼰��Ȩϵ��
		Gauss gau;
		if(gau.Generate_gauss(gaupoi_num)==0) return 0;
		
		//Ԥ���ֲ�ģ�͵�Ч�ն���
		MathMatrix Km(6, 6);
		for(int i=0; i<6; i++)
		{
			double a[6] = {0};
			a[i] = 1.0;
			if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, a, mat, peri_para.horizon_R, peri_para.intrinsic_L, 0)==0) return 0;
			Km.element[0][i] = mat.elas_matrix[0][0];
			Km.element[1][i] = mat.elas_matrix[1][1];
			Km.element[2][i] = mat.elas_matrix[2][2];
			Km.element[3][i] = mat.elas_matrix[0][1];
			Km.element[4][i] = mat.elas_matrix[1][2];
			Km.element[5][i] = mat.elas_matrix[0][2];
		}

		//���¼���
		double Cm[6] = {	known_mat.elas_matrix[0][0],	known_mat.elas_matrix[1][1],	known_mat.elas_matrix[2][2], 
										known_mat.elas_matrix[0][1], known_mat.elas_matrix[1][2], known_mat.elas_matrix[0][2] };

		MathMatrix Dm(6,6);
		Dm = (Km.Transpose()*Km).Inverse()*Km.Transpose();
		for(int i=0; i<6; i++)
		{
			peri_para.acoe[i] = 0.0;
			for(int j=0; j<6; j++) peri_para.acoe[i] += Dm.element[i][j]*Cm[j];
		}

		if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, peri_para.acoe, mat, peri_para.horizon_R, peri_para.intrinsic_L, 1)==0) return 0;

		//----------------------------------------------------------------------------
		//��ȷ�����Ч�ն���
//		for(int i=0; i<6; i++)
//			for(int j=0; j<6; j++)
//			{
//				mat.elas_matrix[i][j] = 0; //����
//				mat.elas_matrix[i][j] += Stiffness_long_range_interaction_Cartesian(i, j, gau.gauss, gau.weight, peri_para.horizon_R, peri_para.intrinsic_L);  //�ѿ����������	
////				mat.elas_matrix[i][j] += Stiffness_long_range_interaction_Spherical(i, j, gau.gauss, gau.weight, peri_para.horizon_R, peri_para.intrinsic_L);  //���������
//			}

		//---------------------------------------------------------------------------
		//�������۹�ʽ����ϵ��������֤һ���ԣ�����Ҫ����Horizon�뾶������Horizon�뾶Ӧ�ô���1������������
//		mat.Compare_coef_by_analysis_formula(stif_nonloc.type, peri_para.horizon_R, peri_para.acoe);
		//---------------------------------------------------------------------------
		//�������۹�ʽ����նȾ��󣬲���֤һ����
//		mat.Compare_matrix_by_analysis_formula(peri_para.horizon_R, peri_para.acoe);

		//---------------------------------------------------------------------------
		//���ݵ��Ծ��󣬷�����ϲ���
		mat.Get_ele_para_by_ela_matrix();

		mats_vec.push_back(mat); //�����������
	}

	//---------------------------------------------------------------------------
	//������в��ϵĸնȾ��󼰶�Ӧ������ģ�������ɱ�
	Print_stiffness_matrix("Material_stiffness_matrix.dat", peri_para.acoe);
	//���ݵ��Ծ��󣬷�����ϲ���
	known_mat.Get_ele_para_by_ela_matrix();
	Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 1);

	return 1;
}
//---------------------------------------------------------------------------
//Generate the material base with damage information
int MatBase::Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para)
{
	MatPro Total_mat(stif_nonloc);
	if(Total_mat.Generate_elas_matrix()==0) return 0; //���ɲ��ϵĸնȾ���

	//---------------------------------------------------------------------------
	MatPro known_mat;
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			known_mat.elas_matrix[i][j] = Total_mat.elas_matrix[i][j]*(1.0-damages.d_crit);
	known_mat.Get_ele_para_by_ela_matrix();

	//---------------------------------------------------------------------------
	//�޸Ĳ���ϵ�������ݷǾֲ�������
	if(peri_para.horizon_R>0.0) 
	{
		MatPro mat(stif_nonloc.type);
		//---------------------------------------------------------------------------
		//���ɸ�������Ԥ���ֲ�ģ�͵�Ч�ն���
		Mesher grid;
		double grid_vol;
		if(grid.Generate_grids_for_effective_stiffness(nonloc_gsize, peri_para.horizon_R, grid_vol)==0) return 0;

		//ȡ��˹�㼰��Ȩϵ��
		Gauss gau;
		if(gau.Generate_gauss(gaupoi_num)==0) return 0;
		
		//Ԥ���ֲ�ģ�͵�Ч�ն���
		MathMatrix Km(6, 6);
		for(int i=0; i<6; i++)
		{
			double a[6] = {0};
			a[i] = 1.0;
			if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, a, mat, peri_para.horizon_R, peri_para.intrinsic_L, 1)==0) return 0;
			Km.element[0][i] = mat.elas_matrix[0][0];
			Km.element[1][i] = mat.elas_matrix[1][1];
			Km.element[2][i] = mat.elas_matrix[2][2];
			Km.element[3][i] = mat.elas_matrix[0][1];
			Km.element[4][i] = mat.elas_matrix[1][2];
			Km.element[5][i] = mat.elas_matrix[0][2];
		}

		//���¼���
		double Cm[6] = {	known_mat.elas_matrix[0][0],	known_mat.elas_matrix[1][1],	known_mat.elas_matrix[2][2], 
										known_mat.elas_matrix[0][1], known_mat.elas_matrix[1][2], known_mat.elas_matrix[0][2] };

		MathMatrix Dm(6,6);
		Dm = (Km.Transpose()*Km).Inverse()*Km.Transpose();
		for(int i=0; i<6; i++)
		{
			peri_para.acoe[i] = 0.0;
			for(int j=0; j<6; j++) peri_para.acoe[i] += Dm.element[i][j]*Cm[j];
		}

		if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, peri_para.acoe, mat, peri_para.horizon_R, peri_para.intrinsic_L, 1)==0) return 0;

		//----------------------------------------------------------------------------
		//��ȷ�����Ч�ն���
//		for(int i=0; i<6; i++)
//			for(int j=0; j<6; j++)
//			{
//				mat.elas_matrix[i][j] = 0; //����
//				mat.elas_matrix[i][j] += Stiffness_long_range_interaction_Cartesian(i, j, gau.gauss, gau.weight, peri_para.horizon_R, peri_para.intrinsic_L);  //�ѿ����������	
////				mat.elas_matrix[i][j] += Stiffness_long_range_interaction_Spherical(i, j, gau.gauss, gau.weight, peri_para.horizon_R, peri_para.intrinsic_L);  //���������
//			}

		//---------------------------------------------------------------------------
		//�������۹�ʽ����ϵ��������֤һ���ԣ�����Ҫ����Horizon�뾶������Horizon�뾶Ӧ�ô���1������������
//		mat.Compare_coef_by_analysis_formula(stif_nonloc.type, peri_para.horizon_R, peri_para.acoe);
		//---------------------------------------------------------------------------
		//�������۹�ʽ����նȾ��󣬲���֤һ����
//		mat.Compare_matrix_by_analysis_formula(peri_para.horizon_R, peri_para.acoe);

		//---------------------------------------------------------------------------
		//���ݵ��Ծ��󣬷�����ϲ���
		mat.Get_ele_para_by_ela_matrix();

		mats_vec.push_back(mat); //�����������
	}

	//---------------------------------------------------------------------------
	//The total siffness matrix (for damage part)
	mats_vec.push_back(Total_mat);

	//---------------------------------------------------------------------------
	//���ݵ��Ծ��󣬷�����ϲ���
	Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 0);
	Print_stiffness_matrix("Material_stiffness_matrix.dat", Total_mat, 1);

	//---------------------------------------------------------------------------
	//A material properties with small damage
	//for(int i=0; i<6; i++)
	//	for(int j=0; j<6; j++)
	//		known_mat.elas_matrix[i][j] = Total_mat.elas_matrix[i][j]*0.80;
	//known_mat.Get_ele_para_by_ela_matrix();
	//mats_vec.push_back(known_mat);
	//Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 1);

	//---------------------------------------------------------------------------
	//������в��ϵĸնȾ��󼰶�Ӧ������ģ�������ɱ�
	Print_stiffness_matrix("Material_stiffness_matrix.dat", peri_para.acoe, 1);

	return 1;
}
//------------------------------------------------------------------------------
//���㳤������Ч�ն�ֵ����ά������ϵ��, ���㹫ʽΪexp(x/L)*x*x
double MatBase::Stiffness_long_range_interaction_Spherical(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len)
{
	double value = 0.0;
	double r, sita, phi;
	for(int i=0; i<(int)weight.size(); i++)
	{
		r = 0.5*(1.0+gauss[i].x)*decayR;
		sita = 0.5*(1.0+gauss[i].y)*PI;
		phi = 0.5*(1.0+gauss[i].z)*2.0*PI;
		value = value + weight[i]*exp(-r/inst_len)*Spherial_coordinates(ni, r, sita, phi)*Spherial_coordinates(nj, r, sita, phi)*r*r*sin(sita); //��������ϵ�»���, ��Ҫ������������:( r, sita, phi )�ֱ��Ӧ( 1, r, rsin(phi) ), ���� r*r*sin(sita)
	}

	value = 0.5*(0.5*decayR*0.5*PI*PI)*value;

	return value;
}
//------------------------------------------------------------------------------
//����������������
double MatBase::Spherial_coordinates(const int &num, const double &r, const double &sita, const double &phi)
{
	double value=0.0;

	switch(num)
	{
	case 0:{ value = r*sin(sita)*cos(phi)*r*sin(sita)*cos(phi); break; }
	case 1:{ value = r*sin(sita)*sin(phi)*r*sin(sita)*sin(phi); break; }
	case 2:{ value = r*cos(sita)*r*cos(sita); break; }
	case 3:{ value = r*sin(sita)*cos(phi)*r*sin(sita)*sin(phi); break; }
	case 4:{ value = r*sin(sita)*sin(phi)*r*cos(sita); break; }
	case 5:{ value = r*sin(sita)*cos(phi)*r*cos(sita); break; }
	default: hout << "����������������ʱ�����������������飡" << endl;
	}

	return value;
}
//------------------------------------------------------------------------------
//���㳤������Ч�ն�ֵ����ά�ѿ�������ϵ��, ���㹫ʽΪe^(x/L)*x*x
double MatBase::Stiffness_long_range_interaction_Cartesian(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len)
{
	double value = 0.0;
	for(int i=0; i<(int)weight.size(); i++)
	{
		double x, y, z;
		x = gauss[i].x*decayR;
		y = gauss[i].y*decayR;
		z = gauss[i].z*decayR;
		double squ_dis = x*x + y*y + z*z; //����֮������ƽ��
		double poi_dis = sqrt(squ_dis);  //����֮��ľ���
		if(poi_dis>decayR+Zero) continue;
		value = value + weight[i]*exp(-poi_dis/inst_len)*Cartesian_coordinates(ni, x, y, z)*Cartesian_coordinates(nj, x, y, z); 
	}

	value = 0.5*decayR*decayR*decayR*value;

	return value;
}
//------------------------------------------------------------------------------
//����ѿ�����������
double MatBase::Cartesian_coordinates(const int &num, const double &x, const double &y, const double &z)
{
	double value=0.0;

	switch(num)
	{
	case 0:{ value = x*x; break; }
	case 1:{ value = y*y; break; }
	case 2:{ value = z*z; break; }
	case 3:{ value = x*y; break; }
	case 4:{ value = y*z; break; }
	case 5:{ value = z*x; break; }
	default: hout << "����ѿ�������ϵ����ʱ�����������������飡" << endl;
	}

	return value;
}
//------------------------------------------------------------------------------
//������в��ϵĸնȾ����ļ���
void MatBase::Print_stiffness_matrix(string print_name, const double acoe[], const int &key)const
{
	ofstream opri;
	if(key==0) opri.open(print_name.c_str());
	else opri.open(print_name.c_str(), ios_base::app);
	//���ϸնȾ������
	opri<<endl;
	opri<<"effective coefficients for the non local model"<<endl;
	opri<<"----------------------------------------------"<<endl;
	for (int i=0; i<6;i++)
	{
		opri<<"coef "<<i<<" = "<<acoe[i]<<endl;
	}
	opri<<"----------------------------------------------"<<endl;
	for(int num=0; num<(int)mats_vec.size(); num++)
	{
		opri << "Material " << num+1 << "  stiffness_matrix:" << endl;
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mats_vec[num].elas_matrix[i][j] << "  ";
			opri << endl;
		}
		opri << endl;

		//���ϵ�����ģ�������ɱȺͼ���ģ�����
		opri << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu23��Nu13 ; Shear Modulus��G12��G23��G13 " << endl;
		opri << mats_vec[num].E11 << " " << mats_vec[num].E22 << " " << mats_vec[num].E33 << "   ";
		opri << mats_vec[num].Nu12 << " " << mats_vec[num].Nu23 << " " << mats_vec[num].Nu13 << "   ";
		opri << mats_vec[num].G12 << " " << mats_vec[num].G23 << " " << mats_vec[num].G13 << "   ";

		opri << endl << endl;
	}

	opri.close();
}
//------------------------------------------------------------------------------
//������ϵĸնȾ����ļ���
void MatBase::Print_stiffness_matrix(string print_name, const double acoe[])const
{
	ofstream opri(print_name.c_str());
	//���ϸնȾ������
	opri<<endl;
	opri<<"effective coefficients for the non local model"<<endl;
	opri<<"----------------------------------------------"<<endl;
	for (int i=0; i<6;i++)
	{
		opri<<"coef "<<i<<" = "<<acoe[i]<<endl;
	}
	opri<<"----------------------------------------------"<<endl;
	for(int num=0; num<(int)mats_vec.size(); num++)
	{
		opri << "Material " << num+1 << "  stiffness_matrix:" << endl;
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mats_vec[num].elas_matrix[i][j] << "  ";
			opri << endl;
		}
		opri << endl;

		//���ϵ�����ģ�������ɱȺͼ���ģ�����
		opri << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu23��Nu13 ; Shear Modulus��G12��G23��G13 " << endl;
		opri << mats_vec[num].E11 << " " << mats_vec[num].E22 << " " << mats_vec[num].E33 << "   ";
		opri << mats_vec[num].Nu12 << " " << mats_vec[num].Nu23 << " " << mats_vec[num].Nu13 << "   ";
		opri << mats_vec[num].G12 << " " << mats_vec[num].G23 << " " << mats_vec[num].G13 << "   ";

		opri << endl << endl;
	}

	opri.close();
}
//------------------------------------------------------------------------------
//������ϵĸնȾ����ļ���
void MatBase::Print_stiffness_matrix(string print_name, MatPro &mat, const int &key)const
{
	ofstream opri;
	if(key==0) opri.open(print_name.c_str());
	else opri.open(print_name.c_str(), ios_base::app);
	//���ϸնȾ������
	opri << "Material  stiffness_matrix:" << endl;
	for(int i=0; i<6; i++)
	{
		for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mat.elas_matrix[i][j] << "  ";
		opri << endl;
	}
	opri << endl;

	//���ϵ�����ģ�������ɱȺͼ���ģ�����
	opri << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu23��Nu13 ; Shear Modulus��G12��G23��G13 " << endl;
	opri << mat.E11 << " " << mat.E22 << " " << mat.E33 << "   ";
	opri << mat.Nu12 << " " << mat.Nu23 << " " << mat.Nu13 << "   ";
	opri << mat.G12 << " " << mat.G23 << " " << mat.G13 << "   ";

	opri << endl << endl;

	opri.close();
}
//------------------------------------------------------------------------------
//Ԥ���ֲ�ģ�͵�Ч�ն���
int MatBase::Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
																const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const
{
	//--------------------------------------------------	
	//��������Ԫ��ż���˹������
	if((int)elements.size()==0) { hout << "Ԥ���ն�����������Ԫ��Ϊ0, ���飡" << endl; return 0; }
	int nume = ((int)elements.size()-1)/2;

	//--------------------------------------------
	//��˹�����
	const int gau_num = (int)weight.size();

	//--------------------------------------------
	//��������Ԫ��˹������
	Point_3D gd_temp(0,0,0);
	vector<Point_3D> gd_main(gau_num, gd_temp);

	for(int count=0; count<gau_num; count++)
	{
		double Nshape[8] = {0};
		
		//�����˹���ֵ���������
		Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

		//�����˹���ֵ�����
		for(int j=0; j<8; j++) 
		{
			gd_main[count].x += Nshape[j]*nodes[elements[nume].nodes_id[j]].x;
			gd_main[count].y += Nshape[j]*nodes[elements[nume].nodes_id[j]].y;
			gd_main[count].z += Nshape[j]*nodes[elements[nume].nodes_id[j]].z;
		}
	}

	//--------------------------------------------------	
	//��������Ч�ն���
	double (*Long_range_stiffness)[6][6] = new double [gau_num][6][6];
	//��������Ч�ն����ʼ��Ϊ��
	for(int count =0; count<gau_num; count++)		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				Long_range_stiffness[count][i][j] = 0.0;

	//--------------------------------------------------	
	//ѭ�����е�Ԫ
	for(int i=0; i<(int)elements.size(); i++)
	{
		//ѭ������Ԫ��˹��������
		for(int count1=0; count1<gau_num; count1++)
		{
			//ѭ���ⵥԪ��˹��������
			for(int count2=0; count2<gau_num; count2++)
			{
				double Nshape[8] = {0};
				//--------------------------------------------
				//�����˹���ֵ���������
				Nshape[0]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[1]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[2]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[3]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[4]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[5]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[6]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[7]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);

				//--------------------------------------------
				//�����˹���ֵ�����
				Point_3D gd_assist(0,0,0);
				for(int j=0; j<8; j++) 
				{
					gd_assist.x += Nshape[j]*nodes[elements[i].nodes_id[j]].x;
					gd_assist.y += Nshape[j]*nodes[elements[i].nodes_id[j]].y;
					gd_assist.z += Nshape[j]*nodes[elements[i].nodes_id[j]].z;
				}				
			
				//--------------------------------------------
				//���㳤������˥������ֵ
				const double x = gd_assist.x-gd_main[count1].x;
				const double y = gd_assist.y-gd_main[count1].y;
				const double z = gd_assist.z-gd_main[count1].z;	

				const double dis_squr = x*x + y*y + z*z;
				const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
				if(poi_dis>decayR||poi_dis<Zero) continue;		//Բ�λ�����
			
				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				double sum = a[0] + a[1]*0.5*(3*cos2sita-1) + a[2]*(2*cos2pha-1)*3*(1-cos2sita) + a[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
										+ a[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + a[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);

				const double gv = exp(-poi_dis/inst_len)*sum*weight[count2]; //˥���Լ�Ȩ�غ���ֵ
//				const double gv = sum*weight[count2]/dis_squr; //˥���Լ�Ȩ�غ���ֵ  (���ڶԱȷ�����ʽ����ֵ��ʽ)

				const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //��Ч������˥������
				const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //˥����������ĶԳ���

				//--------------------------------------------
				//���㳤������Ч�նȾ���
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						Long_range_stiffness[count1][j][k] += temp_gmat[j]*r_comp[k]; //�൱��gv*r_comp[j]*r_comp[k]
			}
		}
	}

	//�����ſɱ�����ʽ��ֵ
	//--------------------------------------------
	//�κ���N�����ĵ��ƫ������
	double diff[3][8];
	diff[0][0]=-0.125;
	diff[0][1]=-diff[0][0];                         
	diff[0][2]=0.125;
	diff[0][3]=-diff[0][2];
	diff[0][4]=-0.125;
	diff[0][5]=-diff[0][4];
	diff[0][6]=0.125;
	diff[0][7]=-diff[0][6];

	diff[1][0]=-0.125;
	diff[1][1]=-0.125;
	diff[1][2]=-diff[1][1];
	diff[1][3]=-diff[1][0];
	diff[1][4]=-0.125;
	diff[1][5]=-0.125;
	diff[1][6]=-diff[1][5];
	diff[1][7]=-diff[1][4];

	diff[2][0]=-0.125;
	diff[2][1]=-0.125;
	diff[2][2]=-0.125;
	diff[2][3]=-0.125;
	diff[2][4]=-diff[2][0];
	diff[2][5]=-diff[2][1];
	diff[2][6]=-diff[2][2];
	diff[2][7]=-diff[2][3];
	//--------------------------------------------------
	//��Ԫ�ڵ��������
	double elenode[8][3];
	for(int j=0; j<8; j++)
	{
		elenode[j][0]=nodes[elements[nume].nodes_id[j]].x;
		elenode[j][1]=nodes[elements[nume].nodes_id[j]].y;
		elenode[j][2]=nodes[elements[nume].nodes_id[j]].z;
	}
	//--------------------------------------------------
	//J����
	double Jmatrix[3][3];
	//������������Ļ�
	for(int j=0; j<3; j++)
		for(int k=0; k<3; k++)
		{
			Jmatrix[j][k]=0;
			for(int m=0; m<8; m++)
			Jmatrix[j][k] += diff[j][m]*elenode[m][k];
		}
	//--------------------------------------------------
	//���J���������ʽ
	double J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
								-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
								+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);

	//--------------------------------------------------
	//�������������С����
	double error_max[6][6], error_min[6][6];	

	//--------------------------------------------------
	//ѭ������Ԫ��˹����������˹���ϵĵ�Ч����������
	const double deco = 0.5*J_val;  //���̵�Ч�����ϵ��
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
		{
			mat.elas_matrix[i][j] = 0;
			for(int count=0; count<gau_num; count++)
			{
				Long_range_stiffness[count][i][j] = Long_range_stiffness[count][i][j]*deco;

				if(count==0) { error_max[i][j] = Long_range_stiffness[count][i][j];  error_min[i][j] = Long_range_stiffness[count][i][j]; }
				else 
				{
					if(error_max[i][j] < Long_range_stiffness[count][i][j]) error_max[i][j] = Long_range_stiffness[count][i][j];
					if(error_min[i][j] > Long_range_stiffness[count][i][j]) error_min[i][j] = Long_range_stiffness[count][i][j];
				}
				
				mat.elas_matrix[i][j] += Long_range_stiffness[count][i][j]*weight[count];
			}
			mat.elas_matrix[i][j] = mat.elas_matrix[i][j]*J_val/grid_vol;
		}

	if(output_key==1)
	{
		//����������
		//for(int count=0; count<gau_num; count++)
		//{
		//	hout << endl << "The " << count << " effective matrix:" << endl;
		//	for(int i=0; i<6; i++)
		//	{
		//		for(int j=0; j<6; j++)
		//			hout << Long_range_stiffness[count][i][j] << "  ";
		//		hout << endl;
		//	}
		//}

		//���������ڼ��
		//hout << endl << "The final error effective matrix:" << endl;
		//for(int i=0; i<6; i++)
		//{
		//	for(int j=0; j<6; j++)
		//		hout << (error_max[i][j]-error_min[i][j])/mat.elas_matrix[i][j] << "  ";
		//	hout << endl;
		//}
	}

	return 1;
}
//------------------------------------------------------------------------------
//Ԥ���ֲ�ģ�͵�Ч�ն���
int MatBase::Estimate_stiffness_long_range_interaction_2D(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
																const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const
{
	//--------------------------------------------------	
	//��������Ԫ��ż���˹������
	if((int)elements.size()==0) { hout << "Ԥ���ն�����������Ԫ��Ϊ0, ���飡" << endl; return 0; }
	int nume = ((int)elements.size()-1)/2;

	//--------------------------------------------
	//��˹�����
	const int gau_num = (int)weight.size();

	//--------------------------------------------
	//��������Ԫ��˹������
	Point_3D gd_temp(0,0,0);
	vector<Point_3D> gd_main(gau_num, gd_temp);

	for(int count=0; count<gau_num; count++)
	{
		double Nshape[4] = {0};
		
		//�����˹���ֵ���������
		Nshape[0]=0.25*(1.0-gauss[count].x)*(1.0-gauss[count].y);
		Nshape[1]=0.25*(1.0+gauss[count].x)*(1.0-gauss[count].y);
		Nshape[2]=0.25*(1.0+gauss[count].x)*(1.0+gauss[count].y);
		Nshape[3]=0.25*(1.0-gauss[count].x)*(1.0+gauss[count].y);

		//�����˹���ֵ�����
		for(int j=0; j<4; j++) 
		{
			gd_main[count].x += Nshape[j]*nodes[elements[nume].nodes_id[j]].x;
			gd_main[count].y += Nshape[j]*nodes[elements[nume].nodes_id[j]].y;
		}
	}

	//--------------------------------------------------	
	//��������Ч�ն���
	double (*Long_range_stiffness)[6][6] = new double [gau_num][6][6];
	//��������Ч�ն����ʼ��Ϊ��
	for(int count =0; count<gau_num; count++)		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				Long_range_stiffness[count][i][j] = 0.0;

	//--------------------------------------------------	
	//ѭ�����е�Ԫ
	for(int i=0; i<(int)elements.size(); i++)
	{
		//ѭ������Ԫ��˹��������
		for(int count1=0; count1<gau_num; count1++)
		{
			//ѭ���ⵥԪ��˹��������
			for(int count2=0; count2<gau_num; count2++)
			{
				double Nshape[4] = {0};
				//--------------------------------------------
				//�����˹���ֵ���������
				Nshape[0]=0.25*(1.0-gauss[count2].x)*(1.0-gauss[count2].y);
				Nshape[1]=0.25*(1.0+gauss[count2].x)*(1.0-gauss[count2].y);
				Nshape[2]=0.25*(1.0+gauss[count2].x)*(1.0+gauss[count2].y);
				Nshape[3]=0.25*(1.0-gauss[count2].x)*(1.0+gauss[count2].y);

				//--------------------------------------------
				//�����˹���ֵ�����
				Point_3D gd_assist(0,0,0);
				for(int j=0; j<4; j++) 
				{
					gd_assist.x += Nshape[j]*nodes[elements[i].nodes_id[j]].x;
					gd_assist.y += Nshape[j]*nodes[elements[i].nodes_id[j]].y;
				}				
			
				//--------------------------------------------
				//���㳤������˥������ֵ
				const double x = gd_assist.x-gd_main[count1].x;
				const double y = gd_assist.y-gd_main[count1].y;

				const double dis_squr = x*x + y*y;
				const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
				if(poi_dis>decayR||poi_dis<Zero) continue;		//Բ�λ�����

				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				double sum = a[0] + a[1]*(2*cos2pha-1) + a[2]*(8*cos2pha*cos2pha-8*cos2pha+1);

//				const double gv = sum*weight[count2]/dis_squr; //constant micro-modulus
				const double gv = exp(-poi_dis/inst_len)*sum*weight[count2]/dis_squr; //˥���Լ�Ȩ�غ���ֵ
//				const double gv = sum*weight[count2]/dis_squr; //˥���Լ�Ȩ�غ���ֵ(���ڶԱȷ�����ʽ����ֵ��ʽ)

				const double r_comp[6] = {x*x, y*y, 0, x*y, y*0, 0*x}; //��Ч������˥������
				const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //˥����������ĶԳ���

				//--------------------------------------------
				//���㳤������Ч�նȾ���
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						Long_range_stiffness[count1][j][k] += temp_gmat[j]*r_comp[k]; //�൱��gv*r_comp[j]*r_comp[k]
			}
		}
	}

	//�����ſɱ�����ʽ��ֵ
	//--------------------------------------------
	//�κ���N�����ĵ��ƫ������
	double diff[2][4];
	diff[0][0]=-0.25;
	diff[0][1]=-diff[0][0];                         
	diff[0][2]=0.25;
	diff[0][3]=-diff[0][2];

	diff[1][0]=-0.25;
	diff[1][1]=-0.25;
	diff[1][2]=-diff[1][1];
	diff[1][3]=-diff[1][0];

	//--------------------------------------------------
	//��Ԫ�ڵ��������
	double elenode[4][2];
	for(int j=0; j<4; j++)
	{
		elenode[j][0]=nodes[elements[nume].nodes_id[j]].x;
		elenode[j][1]=nodes[elements[nume].nodes_id[j]].y;
	}
	//--------------------------------------------------
	//J����
	double Jmatrix[2][2];
	//������������Ļ�
	for(int j=0; j<2; j++)
		for(int k=0; k<2; k++)
		{
			Jmatrix[j][k]=0;
			for(int m=0; m<4; m++)
			Jmatrix[j][k] += diff[j][m]*elenode[m][k];
		}
	//--------------------------------------------------
	//���J���������ʽ
	double J_val = Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0];

	//--------------------------------------------------
	//�������������С����
	double error_max[6][6], error_min[6][6];	

	//--------------------------------------------------
	//ѭ������Ԫ��˹����������˹���ϵĵ�Ч����������
	const double deco = 0.5*J_val;  //���̵�Ч�����ϵ��

	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
		{
			mat.elas_matrix[i][j] = 0;
			for(int count=0; count<gau_num; count++)
			{
				Long_range_stiffness[count][i][j] = Long_range_stiffness[count][i][j]*deco;

				if(count==0) { error_max[i][j] = Long_range_stiffness[count][i][j];  error_min[i][j] = Long_range_stiffness[count][i][j]; }
				else 
				{
					if(error_max[i][j] < Long_range_stiffness[count][i][j]) error_max[i][j] = Long_range_stiffness[count][i][j];
					if(error_min[i][j] > Long_range_stiffness[count][i][j]) error_min[i][j] = Long_range_stiffness[count][i][j];
				}
				
				mat.elas_matrix[i][j] += Long_range_stiffness[count][i][j]*weight[count];
			}
			mat.elas_matrix[i][j] = mat.elas_matrix[i][j]*J_val/grid_vol;
		}

	if(output_key==1)
	{
		//����������
		//for(int count=0; count<gau_num; count++)
		//{
		//	hout << endl << "The " << count << " effective matrix:" << endl;
		//	for(int i=0; i<6; i++)
		//	{
		//		for(int j=0; j<6; j++)
		//			hout << Long_range_stiffness[count][i][j] << "  ";
		//		hout << endl;
		//	}
		//}

		//���������ڼ��
		//hout << endl << "The final error effective matrix:" << endl;
		//for(int i=0; i<6; i++)
		//{
		//	for(int j=0; j<6; j++)
		//		hout << (error_max[i][j]-error_min[i][j])/mat.elas_matrix[i][j] << "  ";
		//	hout << endl;
		//}
	}

	return 1;
}
//===========================================================================
