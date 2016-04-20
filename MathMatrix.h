//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MathMatrix.h
//OBJECTIVE:	Mathematical matrix operations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef MATHMATRIX_H
#define MATHMATRIX_H

#include<cmath>
#include<iomanip>
#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

const double Zero = 1.0E-8;

//---------------------------------------------------------------------------
//������
class MathMatrix
{
public:
	//-----------------------------------------------
	//���ݳ�Ա
	vector<vector<double> > element;

	//-----------------------------------------------
	//��Ա����
	//���캯��
	MathMatrix(){};
	MathMatrix(const int&, const int&);							//�����СΪM*N��ֵΪ��ľ���
	MathMatrix(const vector<double>&);						//һά��������
	MathMatrix(const vector<vector<double> >&);			//��ά��������					
	MathMatrix(const double *, const int&);					//һά���鹹��
	MathMatrix(const double *, const int&, const int&);	//��ά���鹹�죬ʵ��Ϊ��&vec[0][0],m,n��
	//-----------------------------------------------
	void Evalu(const double *, const int&);						//һά���鸳������ 
	void Evalu(const double *, const int& , const int&);	//��ά���鸳������ʵ��Ϊ��&vec[0][0],m,n�� 
	void operator=(const vector<double>&);					//���ظ�ֵ��������һά������������
	void operator=(const vector<vector<double> >&);	//���ظ�ֵ����������ά������������
	MathMatrix operator*(const MathMatrix&);				//���س˷�������������˾���
	MathMatrix operator*(const double&);						//���س˷��������������ʵ����
	MathMatrix operator+(const MathMatrix&);				//���ؼӷ�������������Ӿ���			
	double operator+(const double&);								//���ؼӷ���������1*1�����ʵ����			
	MathMatrix operator-(const MathMatrix&);				//���ؼ��������������������
	MathMatrix Inverse(void);					//��������
	MathMatrix Transpose(void);			//����ת��
	int RowN(void);								//���������
	int CalN(void);									//���������
	int Symmetry(void);							//�ж�����Գ�
	//-----------------------------------------------------------------------------
	MathMatrix GetCol(int cn);                  //�������һ�У�����תΪ����

	//���������������
	friend ostream& operator<<(ostream& o, const MathMatrix& matrix);
};
//���������������
inline ostream& operator<<(ostream& o,const MathMatrix& matrix)
{
	for (int i=0; i<int(matrix.element.size()); i++)
	{
		for(int j=0; j<int(matrix.element[i].size()); j++)
		{
//			if (fabs(matrix.element[i][j])<1.0E-8)
//			{
//				o <<setw(8) << setfill(' ') << 0 << " ";
//			}
//			else
//			{
				o <<setw(16) << setfill(' ') << matrix.element[i][j] << " ";
//			}
		}
		o << endl;
	}
	return o;
}
//---------------------------------------------------------------------------
#endif
//===========================================================================
