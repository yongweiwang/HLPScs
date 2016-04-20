//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Gauss.h
//OBJECTIVE:	Compute the gaussian points
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GAUSS_H
#define GAUSS_H

#include<iostream>
#include<vector>
#include<cmath>
#include "Hns.h"
#include "Fem_3D.h"
using namespace hns;

const double PI = 3.1415926535897932384626433832795;

//----------------------------------------------------------

class Gauss
{
	public:
		//���ݱ���
		vector<Node> gauss; 
		vector<double> weight;

		//���캯����
		Gauss(){};
	
		//��Ա����
		//������ά��˹�ڵ�
		int Generate_gauss(const int &ipre);
		//���ɶ�ά��˹�ڵ�
		int Generate_gauss_2D(const int &ipre);

	private:

		//���ݱ�����
		int precision;

		//�������ɸ�˹������
		int Generate_gauss_array(vector<double> &gauss, vector<double> &weight)const;
		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};//-----------------------------------------------------
#endif
//===========================================================================
