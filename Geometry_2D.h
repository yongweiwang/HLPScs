//====================================================================================
//SOFTWARE:	Non-local Elasto-plastic Continuum Analysis (NECA)
//CODE FILE:	Geometry_2D.h
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
using namespace hns;

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//�����ά�ռ����
class Point_2D
{
	public:

		//���ݳ�Ա
		double x, y;
		int flag;

		//���캯��
		Point_2D(){};
		Point_2D( double px, double py );

		//��Ա����
        Point_2D operator+( Point_2D &pt );
        Point_2D operator+( const Point_2D &pt );
        Point_2D operator+( double d );
        Point_2D operator-( double d );
        Point_2D operator-( Point_2D &pt );
        Point_2D operator-( const Point_2D &pt );
        Point_2D operator*( double d );
        Point_2D operator/( double d );
		bool operator==( Point_2D &pt );
		bool operator!=( Point_2D &pt );
        double distance_to(const Point_2D &pt)const;
		double distance_to(const double &px, const double &py)const;
};

//---------------------------------------------------------------------------
//�����ά��(��)��
class Line_2D
{
	public:

		//���ݳ�Ա
		 Point_2D point[2];	//��ʾ�߶ε������˵�����
		 double len;					//��ʾ�߶εĳ���
		 bool virtual_line;			//��ʾ�߶��Ƿ��϶�Ϊ���߶�(fase:���߶Σ�����Ϊһ�㣩; true:�����߶�)


		//���캯��
		Line_2D(){};
		Line_2D(Point_2D p0, Point_2D p1);
		
		//��Ա����
		double length();	//�߶εĳ���
        double distance_point_to_line(const Point_2D *point_temp)const;    //�㵽ֱ�ߵľ���
        double distance_point_to_line(const Point_2D &point_temp)const;    //�㵽ֱ�ߵľ���
        double distance_point_to_line(double dx, double dy)const;    //�㵽ֱ�ߵľ���
        double path_point_to_line(const Point_2D &point_temp)const;    //�㵽ֱ�ߵ�·�����жϵ�����ֱ�ߵ�һ�໹����һ��
		int contain(const Point_2D &point_temp)const;    //�ж��߶ΰ���һ����
		int contain(const double &px, const double &py)const;    //�ж��߶ΰ���һ����

	private:
				
		//���ݳ�Ա
        double A, B, D ;       //ֱ�߷��̵�ϵ�� Ax+By+D=0
};

//---------------------------------------------------------------------------
//���������
class Rectangle
{
	public:

		//���ݳ�Ա
		Point_2D point[4];	//���������ĸ��ǵ㣬�����½ǿ�ʼ˳ʱ��ת
		double length;			//�������ĳ�
		double width;				//�������Ŀ�
		double area;				//��ʾ���ε����
		int virtual_rect;			//��ʾ�����Ƿ��϶�Ϊ�����(0:�����; 1:�������)

		//���캯��
		Rectangle(){};
		Rectangle (Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3);
		Rectangle (Point_2D p0, double len, double wid);
		Rectangle (Point_2D p0, Point_2D p2);

		//��Ա����
		double calculate_area();			//������ε����
		int contain(const Point_2D &poi)const;		//�жϾ��������Ƿ������һ���㣨����ֵ��0û������1������
		int contain_in(const Point_2D &poi)const;//�ж�һ�����Ƿ񱻰����ھ��������ڲ�(�߽����)������ֵ��0û������1������
		int contain(const Line_2D &line)const;			//�жϾ��������Ƿ������һ���߶Σ�����ֵ��0û������1������
		int contain(const Rectangle &rect)const;		//�жϾ��������Ƿ��������һ���Σ�����ֵ��0û������1������
		int contain_in(const Line_2D &line)const;	//�ж�һ���߶��Ƿ񱻰����ھ����ڲ����ھ��б߽��ϲ��㣨����ֵ��0û������1������ 
		int overlap(const Line_2D &line)const;			//�ж��߶��Ƿ�;����ཻ���ཻ�����߶ε������˵㶼�ھ��εı߽��ϣ�������һ���˵��ھ�����
																					//���������һ���˵��ھ��α߽��ϣ���һ�˵��ھ����ⲿ�����ཻ)(����ֵ��0���ཻ��1�ཻ)
		void output_parameters();		//������ε�ȫ����Ϣ�������ĸ��������꣬���������������α�ǣ�
		int rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2);	 //�������Rect1��Rect2���ص�������
		int make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids);	//���������ڲ�С����Ϊ�������Ź���ָ�
};
//---------------------------------------------------------------------------
//��¼�������Ԫ������1��2������ཻ��Ԫ���ཻС����
struct Relative_Ele_Rect
{	
	vector<int> relative_ele_num[2];
	vector<Rectangle> relative_rect[2];
};
#endif
//===========================================================================
