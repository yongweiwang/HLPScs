//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatBase.h
//OBJECTIVE:	Generate material database
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================
#ifndef MATBASE_H
#define MATBASE_H

#include "Input_Reader.h"
#include "MatPro.h"
#include "Gauss.h"
#include "Mesher.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class MatBase
{
	public:

		//���ݳ�Ա
		vector<MatPro> mats_vec;

		//���캯��
		MatBase(){};

		//��Ա����
		//Generate the material base
		int Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, struct Peri_para &peri_para);
		//Generate the material base with damage information 
		int Generate_matbase(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para);
		//Generate the material base with damage information for 2D
		int Generate_matbase_2D(const struct Stif_nonloc &stif_nonloc, const struct Grid_size &nonloc_gsize, const int &gaupoi_num, const struct Damage &damages, struct Peri_para &peri_para);

	private:
		//��Ա����
		//���㳤������Ч�ն�ֵ����ά������ϵ��
		double Stiffness_long_range_interaction_Spherical(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len);
		//���㳤������Ч�ն�ֵ����ά�ѿ�������ϵ��
		double Stiffness_long_range_interaction_Cartesian(const int &ni, const int &nj, const vector<Node> &gauss, const vector<double> &weight, const double &decayR, const double &inst_len);
		//����������������
		double Spherial_coordinates(const int &num, const double &r, const double &sita, const double &phi);
		//����ѿ�����������
		double Cartesian_coordinates(const int &num, const double &x, const double &y, const double &z);
		//������в��ϵĸնȾ���
		void Print_stiffness_matrix(string print_name, const double acoe[])const;
		void Print_stiffness_matrix(string print_name, const double acoe[], const int &key)const;
		void Print_stiffness_matrix(string print_name, MatPro &mat, const int &key)const;
		//Ԥ���ֲ�ģ�͵�Ч�ն���
		int Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
													const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const;
		int Estimate_stiffness_long_range_interaction_2D(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, const vector<Node> &gauss, 
													const vector<double> &weight, const double a[], MatPro &mat, const double &decayR, const double &inst_len, const int &output_key)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
