//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MatPro.h
//OBJECTIVE:	Material properties
//AUTHOR:		Fei Han; Yan Azdoud
//E-MAIL:			fei.han@kaust.edu.sa;  yan.azdoud@kaust.edu.sa
//====================================================================================

#ifndef MATPRO_H
#define MATPRO_H

#include "Input_Reader.h"
#include "Gauss.h"
#include "MathMatrix.h"
#include "Hns.h"

using namespace hns;

//---------------------------------------------------------------------------
class MatPro
{
	public:	
	
		//���ݳ�Ա
		string type_val;   //��ά�����������ͣ�Ŀǰ��������("Isotropic"������ͬ��; "Transverse"����۸���ͬ��; "Orthotropic"��������������)
		double E11, E22, E33, Nu12, Nu23, Nu13, G12, G23, G13;	//���ϵĵ���ģ�������ɱȣ�����ģ��
		double elas_matrix[6][6];		//Dij ���Ծ���

		//���캯����
		MatPro(string itype="Isotropic"){ type_val = itype; }
		MatPro(const struct Stif_nonloc &stif_nonloc);

		//---------------------------------------------
		//��Ա������
		//���ø���ͬ�Բ��ϲ����ĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double &E);
		//���ú�۸���ͬ�Բ��ϵĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double &E1, const double &E2, const double &Nu2);
		//���������������Բ��ϵĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double &iE1, const double &iE2, const double &iE3, const double &iNu12, const double &iNu23, const double &iNu13);
		//����Dij���Ծ��󣬴˵��Ծ����Ӧ��Ӧ������Ϊ��e11,e22,e33,2*e12,2*e23,2*e31��ת��
		int Generate_elas_matrix();
		int Generate_elas_matrix_2D();
		//����Dij���Ծ��󣬷�����ϲ���
		int Get_ele_para_by_ela_matrix();
		//�������۹�ʽ����ϵ��������֤һ���ԣ�����Ҫ����Horizon�뾶������Horizon�뾶Ӧ�ô���1������������
		void Compare_coef_by_analysis_formula(const int mat_type, const double &decayR, const double acoe[])const;
		//�������۹�ʽ����նȾ��󣬲���֤һ����
		void Compare_matrix_by_analysis_formula(const double &decayR, const double decay_acoe[])const;

		void print()const;
};
//------------------------------------------------
#endif
//===========================================================================
