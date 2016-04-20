//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	MathMatrix.cpp
//OBJECTIVE:	Mathematical matrix operations
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//���캯��

//���캯����M*N����ʼ����
MathMatrix::MathMatrix(const int& M, const int& N)
{
	vector<double> temp(N,0.0);
	element.assign(M,temp);
}
MathMatrix::MathMatrix(const vector<double>& vec2D)
{
	int N = int(vec2D.size());
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec2D[i];
	}
}
MathMatrix::MathMatrix(const vector<vector<double> >& vec2D)
{
	element = vec2D;
}
//���캯�������鹹��1*N��
MathMatrix::MathMatrix(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "��������д�N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//���캯�������鹹��M*N��
MathMatrix::MathMatrix(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "��������д�M<=0��N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//һά���鸳������ 
void MathMatrix::Evalu(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "��������д�N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//---------------------------------------------------------------------------
//��ά���鸳������ 
void MathMatrix::Evalu(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "��������д�M<=0��N<=0 " << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//���ظ�ֵ��������һά������������ 
void MathMatrix::operator=(const vector<double>& vec)
{
	element.assign(1,vec);
}
//���ظ�ֵ����������ά������������
void MathMatrix::operator=(const vector<vector<double> >& vec)
{
	element = vec;
}
//---------------------------------------------------------------------------
//���س˷�������������˾���
MathMatrix MathMatrix::operator*(const MathMatrix& matrix)
{
	int L = 0;
	if (element[0].size() != matrix.element.size())
	{
		cout << "��������ˣ��кŲ������кţ�" << endl;
		exit(0);
	}
	else
	{
		L = int(matrix.element.size());
	}
	//�������M�к���N�к�
	int M = int(element.size());
	int N = int(matrix.element[0].size());
	//����������Tem_mat
	MathMatrix Tem_mat(M,N);
	//�ۼ���
	double conter;

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			conter = 0.0;
			for(int k=0; k<L; k++)
			{
				conter = conter + element[i][k] * matrix.element[k][j];
			}
			Tem_mat.element[i][j] = conter;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//���س˷��������������ʵ����
MathMatrix MathMatrix::operator*(const double& R)
{
	//�������M�к���N�к�
	int M = int(element.size());
	int N = int(element[0].size());
	//����������Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] * R;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//���ؼӷ�������������Ӿ���
MathMatrix MathMatrix::operator+(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "������������кŻ��кŲ���ȣ�" << endl;
		exit(0);
	}
	//�������M�к���N�к�
	int M = int(element.size());
	int N = int(element[0].size());
	//����������Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] + matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//���ؼӷ��������������ʵ����
double MathMatrix::operator+(const double& R)
{
	//�������M�к���N�к�
	int M = int(element.size());
	int N = int(element[0].size());
	if (M!=N||M==1)
	{
		cout << "������1*1�ģ����ܺ�һ��ʵ����ӣ�" << endl;
		exit(0);
	}
	//����������
	double R2 = element[0][0] + R;

	return R2;
}
//---------------------------------------------------------------------------
//���ؼ��������������������
MathMatrix MathMatrix::operator-(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "������������кŻ��кŲ���ȣ�" << endl;
		exit(0);
	}
	//�������M�к���N�к�
	int M = int(element.size());
	int N = int(element[0].size());
	//����������Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] - matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//�ж��Գƾ��� ���-1�����Ƿ������0�����Գƣ����1���Գ�
int MathMatrix::Symmetry(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	if (M==N)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=i+1; j<M; j++)
			{
				if (fabs(element[i][j]-element[j][i])>Zero)
				{
					return 0;
				}
			}
		}
	}
	else
	{
		return -1;
	}
	return 1;
}
//---------------------------------------------------------------------------
//��������
MathMatrix MathMatrix::Inverse(void)
{
	vector<vector<double> > vec = element;
	//--------------------------------------------------
	//�ж��Գ���������

	//--------------------------------------------------
	//�������
	int stRank = int(vec.size());
	vector<double> b(stRank,0.0);
	for(int k=0; k<stRank; k++)
	{
		double w= vec[0][0];
		int m = stRank - k -1;
		for(int i=1; i<stRank; i++)
		{
			double g = vec[i][0];
			b[i] = g / w;
			if (i<=m)
			{
				b[i] = -b[i];
			}
			for(int j=1; j<=i; j++)
			{
				vec[i-1][j-1] = vec[i][j] + g * b[j];
			}
		}
		vec[stRank-1][stRank-1] = 1.0 / w;
		for(int i= 1; i<stRank; i++)
		{
			vec[stRank-1][i-1] =  b[i];
		}
	}
	for(int i=0; i<stRank-1; i++)
	{
		for(int j = i+1; j<stRank; j++)
		{
			vec[i][j] = vec[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//����ת��
MathMatrix MathMatrix::Transpose(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	vector<vector<double> > vec;
	vector<double> temp(M);
	vec.assign(N,temp);

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			vec[i][j] = element[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//�������һ�У�����תΪ����
MathMatrix MathMatrix::GetCol(int cn){
	MathMatrix matrix(int(element.size()),1);
	for(int i=0;i<int(element.size());i++)
		matrix.element[i][0]=element[i][cn];
	return matrix;
}                  
//---------------------------------------------------------------------------
//���������
int MathMatrix::RowN(void)
{
	return int(element.size());
}
//---------------------------------------------------------------------------
//���������
int MathMatrix::CalN(void)
{
	return int(element[0].size());
}
//============================================================================
