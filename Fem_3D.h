//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Fem_3D.h
//OBJECTIVE:	Basic finite element computation classes, including notes and elements, etc.
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef FEM_3D_H
#define FEM_3D_H

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;


//----------------------------------------------------------------------------
// class containing the individual elementary matrices
class Matrices
{
public:

	double Matelem[24][24];
	vector<double> Relative_mat;
};

//---------------------------------------------------------------------------
//����ڵ���
class Node
{	
	public:
		double x, y, z; //�ڵ������ֵ

		int type;	//0�ڵ㣬1�߽���㣬2�߽��ߵ㣬3�ǵ�
		bool DG_node;						//the node for Discontinuous Galerkin computation
		vector<int> relative_nods;	//��ؽڵ��������ڵ��ţ�
		vector<int> relative_eles;	//��ص�Ԫ��������Ԫ��ţ�

		//���캯��
		Node(const double ix=0, const double iy=0, const double iz=0);
		Node(const double ix, const double iy, const double iz, const int itype);
		
		//��Ա������
		double	distance_to(const Node& n); //����
		double	distance_to(Node& n); //����

};
//---------------------------------------------------------------------------
//��Ԫ��ͷ�ļ���
class Element
{
	public:	
		//Data Member
		int type;	//��ʾ��Ԫ����״����(��������xyz�� x ά�ȣ�y��Ԫ�Ľڵ������z��Ԫ���κ����ݴ�)��
						//���磬121: һά���ڵ�(�߶�)�����κ�����231: ��ά���ڵ�(������)�����κ�����241: ��ά�Ľڵ�(�ı���)�����κ�����
						//341: ��ά�Ľڵ�(������)�����κ�����361: ��ά���ڵ�(������)�����κ�����381����ά�˽ڵ�(������)�����κ���
		int mat;	//��ʾ��Ԫ�Ĳ������ԣ�
		int flag;   //��������Ԫ�仯�ı�ǣ�
		bool DG_elem;						//The element for Discontinuous Galerkin computation
		double Damage_energy;		//The damaged energy of this element for continuum damage mechanics
		double Dissipative_energy;	//The dissipative energy of this element for bond broken in peridynamics
		double dam;							//The average damaged value in the center of an element
		vector<int> nodes_id;
		vector<int> relative_eles;	//��ص�Ԫ��������Ԫ��ţ�
		
		//Constructor
		Element(){};
};
//---------------------------------------------------------------------------
//���������嵥Ԫ
class Hexahedron
{
	public:
        int nodesId[8];	 //eight nodes;
        int materialId;
        int facesId[6];

		//���캯����
		Hexahedron(){};
        Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
		
		//member function
		//Calculate the volume of a Hexahedron
		double cal_hex_volume(const vector<Node> &nodes)const;
	private:
		//����3���ڵ�������ɵ�3��3����ʽ��ֵ
		double cal_det(const Node &n1, const Node &n2, const Node &n3)const ;
		//����4�ڵ������嵥Ԫ�����(�ڵ�˳���������ֶ���
		double cal_tet_volume(const Node &node1, const Node &node2, const Node &node3, const Node &node4)const;
};
//---------------------------------------------------------------------------
//������ͷ�ļ���
class Side
{
public:	
	int type;	//��ʾ���ߵ����ͣ�0�����߶Σ�1�������ߣ�2�����ɱ��ߣ�3���̶�����
	vector<int> nodes_id;
};
#endif
//===========================================================================
